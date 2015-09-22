var config = require('./config.json');
var http = require('http');
var https = require('https');
var _ = require('lodash');

function phacts_url( cmd, key, iri ){
  var retval = config['phacts_host'] + config['phacts_base'] + cmd;
  retval += '?app_id=' + config['app_id'] + '&app_key=' + config['app_key'];
  retval += '&' + key + '=' + encodeURIComponent( iri );
  retval += '&_format=json&_pageSize=all';

  return retval;
}

http.createServer(function (req, res) {
  if ( !_.startsWith(req.url, '/drug/') && !_.startsWith(req.url, '/protein/') )
  {
    res.writeHead(404, {'Content-Type': 'text/plain'});
    res.end('404 Not Found');
    return;
  }

  if ( req.url === '/drug/' || req.url === '/protein/' )
  {
    res.writeHead(400, {'Content-Type': 'text/plain'});
    res.end('400 Bad Request');
    return;
  }

  console.log( "Got request: " + req.url );

  if ( _.startsWith(req.url, '/drug/') )
    handle_drug_request( req, res, req.url.substring(6) );
  else
    handle_protein_request( req, res, req.url.substring(9) );
}).listen(config.port);

function handle_protein_request( req, res, iri ){
  MapURI( iri, 'http://www.conceptwiki.org/', function(mapped){
    if ( mapped === null ){
      send200(res, []);
      return;
    }

    if ( mapped.hasOwnProperty('error') ){
      handle_openphacts_error( res, mapped['error'] );
      return;
    }

    var url = phacts_url('target/pharmacology/pages', 'uri', mapped );
    download( url, function(data){
      if ( data.hasOwnProperty('error') ){
        handle_openphacts_error( res, data['error'] );
        return;
      }

      if ( !obj_has_structure( data, ["result", "items"] ) ){
        handle_openphacts_error( res, "OpenPHACTS output did not have the expected structure" );
        return;
      }

      var items = data['result']['items'];

      if ( config['debug'] )
        items = _.chunk(items,25)[0];

      res.targetcnt = items.length;
      res.sparqlcnt = 0;
      res.triples = [];

      for ( var i = 0; i < items.length; i++ ){
        var item = items[i];

        if ( !obj_has_structure( item, ["hasMolecule", "_about"] ) ){
          res.targetcnt--;
          continue;
        }

        res.chemblcnt++;
        MapURI( item['hasMolecule']['_about'], 'http://purl.obolibrary.org/obo/', function(chebi){
          if ( chebi !== null ){
            sparql_query( res, chebi );
          }
          else {
            res.sparqlcnt++;
            if ( res.sparqlcnt >= res.targetcnt ){
              got_all_sparqls( res );
            }
          }
        });
      }

      if ( res.sparqlcnt >= res.targetcnt ){
        got_all_sparqls( res );
      }
    });
  });
}

function handle_drug_request( req, res, id ){
  if ( _.startsWith(id, 'CHEBI') )
    id = "http://purl.obolibrary.org/obo/" + id;

  var PharmURL = phacts_url( 'compound/pharmacology/pages', 'uri', id );

  download( PharmURL, function(data){
    if ( data.hasOwnProperty('error') )
    {
      handle_openphacts_error( res, data['error'] );
      return;
    }

    try{
      var items = data['result']['items'];
    }
    catch(e){
      handle_openphacts_error( res, "OpenPHACTS result did not have the expected format" );
      return;
    }

    var chembls = [];

    for ( var i = 0; i < items.length; i++ ){
      if ( !obj_has_structure( items[i], ['hasAssay', 'hasTarget', '_about'] ) )
        continue;

      chembls.push(items[i]['hasAssay']['hasTarget']['_about']);
    }

    chembls = _.uniq(chembls);

    if ( config['debug'] )
      chembls = _.chunk(chembls,5)[0];

    res.triples = [];
    res.sparqlcnt = 0;
    res.targetcnt = chembls.length;

    for ( var i = 0; i < chembls.length; i++ ){
      MapURI( chembls[i], 'http://www.uniprot.org/', function(mapped){
        sparql_query( res, mapped );
      });
    }

    console.log( "Got " + chembls.length + " non-uniprot IDs.  Converting to uniprot..." );
  });
}

function sparql_query( res, o ){
  if ( o === null || o.hasOwnProperty('error') )
  {
    res.sparqlcnt++;

    if ( res.sparqlcnt === res.targetcnt )
      got_all_sparqls( res );

    return;
  }

  var url = config['sparql'] + encodeURIComponent('SELECT ?s ?p ?o WHERE {?s ?p ?o . FILTER(?o=<'+o+'>)}');

  download(url, function(trips){
    res.triples.push(trips);
    res.sparqlcnt++;

    console.log( 'We have made '+res.sparqlcnt+' (of '+res.targetcnt+') SPARQL queries...' );

    if ( res.sparqlcnt >= res.targetcnt )
      got_all_sparqls( res );
  });
}

function got_all_sparqls( res ){
  var trips = [];

  for ( var i = 0; i < res.triples.length; i++ ){
    if ( res.triples[i].hasOwnProperty('error') ){
      handle_openphacts_error( res, res.triples[i]['error'] );
      return;
    }

    if ( !obj_has_structure( res.triples[i], ["results", "bindings"] ) ){
      handle_openphacts_error( res, 'One of the triple-store response-sets had unexpected structure');
      return;
    }

    trips = _.union( trips, res.triples[i]['results']['bindings'] );
  }

  var unique_trips = [];
  for ( var i = 0; i < trips.length; i++ ){
    if ( !triple_present( trips[i], unique_trips ) ){
      unique_trips.push(trips[i]);
    }
  }

  send200( res, unique_trips );
}

function triple_present( triple, triples ){
  var cmp_trips = function(x,y){
    return x['s']['value']===y['s']['value'] && x['p']['value']===y['p']['value'] && x['o']['value']===y['o']['value'];
  }

  for( var i = 0; i < triples.length; i++ ){
    if ( cmp_trips(triple, triples[i]) )
      return true;
  }

  return false;
}

function MapURI( id, target_base, callback ){
  var url = phacts_url( 'mapUri', 'Uri', id );

  download( url, function(res){
    if ( res.hasOwnProperty('error') ){
      callback(res);
      return;
    }

    if ( !obj_has_structure( res, ['result', 'primaryTopic', 'exactMatch'] ) ){
      callback({'error': "OpenPHACTS MapURI sent a result with unexpected format"});
      return;
    }

    var items = res['result']['primaryTopic']['exactMatch'];

    for ( var i = 0; i < items.length; i++ ){
      if ( _.startsWith( items[i], target_base ) ){
        callback( items[i] );
        return;
      }
    }

    callback( null );
  });
}

function obj_has_structure( o, arr ){
  var optr = o;

  for ( var i = 0; i < arr.length; i++ ){
    if ( !optr.hasOwnProperty(arr[i]) )
      return false;

    optr = optr[arr[i]];
  }

  return true;
}

function send200( res, what ){
console.log("Debug:send200");
console.log("------");
console.log(JSON.stringify(what,null,2));
console.log("------");
  res.writeHead(200, {'Content-Type': 'text/plain'});
  res.end(JSON.stringify(what, null, 2));
}

function handle_openphacts_error( res, e ){
  res.writeHead(500, {'Content-Type': 'text/plain'});
  res.end('500 Internal Server Error\n\nError communicating with OpenPHACTS\nDetails:\n' + e);
}

function download(url, callback) {
  console.log( "Connecting to " + url );
  var handler;

  if ( _.startsWith(url, "https") )
    handler = https;
  else
    handler = http;

  handler.get(url, function(res) {
    var data = "";
    res.on('data', function (chunk) {
      data += chunk;
    });
    res.on("end", function() {
      try{
        callback(JSON.parse(data));
      }
      catch(e){
        callback({'error': ""+e+"\n\nOpenPHACTS sent a non-JSON response:\n"+data});
      }
    });
  }).on("error", function(e) {
    callback({'error': e});
  });
}
