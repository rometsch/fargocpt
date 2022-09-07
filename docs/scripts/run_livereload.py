import socket
import errno
from livereload import Server, shell

import os

def main():
    print(os.getcwd())
    server = Server()
    server.watch('source/**.rst', shell('make html'), delay=1)
    server.watch('source/*.md', shell('make html'), delay=1)
    server.watch('source/*.md', shell('make html'), delay=1)
    for root, dirs, _ in os.walk("source"):
        for d in dirs:
            server.watch(f"{root}/{d}/*.md", shell('make html'), delay=1)
    # server.watch('source/Numerics/*.md', shell('make html'), delay=1)
    server.watch('source/**.py', shell('make html'), delay=1)
    server.watch('_static/**', shell('make html'), delay=1)
    server.watch('_templates/**', shell('make html'), delay=1)

    port = get_next_free_port(5500)
    server.serve(root='build/html', port=port)


def is_port_free(port):
    """ Check whether the port is free to use on localhost.

    Parameters
    ----------
    port: int
        Port number to check.

    Returns
    -------
    bool:
        True if free, false if already in use.
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    try:
        s.bind(("127.0.0.1", port))
        rv = True
    except socket.error as e:
        if e.errno == errno.EADDRINUSE:
            rv = False
        else:
            # something else raised the socket.error exception
            raise

    s.close()

    return rv


def get_next_free_port(port):
    """ Find the next free port on localhost.
    
    Recursively test all port numbers with increments of 1 
    starting from number 'port'.
    
    Parameters
    ----------
    port: int
        Start looking from this port number upwards.
        
    Returns
    -------
    int
        Next free port number.
    """
    if is_port_free(port):
        return port
    else:
        return get_next_free_port(port + 1)
if __name__ == '__main__':
    main()