
import sys, traceback
import BaseHTTPServer
import CGIHTTPServer

protocol="HTTP/1.0"
HandlerClass = CGIHTTPServer.CGIHTTPRequestHandler
ServerClass = BaseHTTPServer.HTTPServer


port = 50005
server_address = ('127.0.0.1', port)
HandlerClass.protocol_version = protocol

try:
    httpd = ServerClass(server_address, HandlerClass)
except:
    #traceback.print_exc()
    print '*** Error starting server, port %s may be in use'  % port
    sys.exit(1)

try:
    sa = httpd.socket.getsockname()
    print "\nServing HTTP on", sa[0], "port", sa[1], "..."
    print 'Use Ctrl-C to shut down server'
    print ' '
    print 'Point your browser to http://localhost:%s/' \
          % port
    print ' '
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print "Server shut down"

except:
    traceback.print_exc()
    print '*** Error starting server'
    sys.exit(1)

