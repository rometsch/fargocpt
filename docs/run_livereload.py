#!.venv/bin/python
from livereload import Server, shell

def main():
    server = Server()
    server.watch('source/*.rst', shell('make html'), delay=1)
    server.watch('source/*.md', shell('make html'), delay=1)
    server.watch('source/*.py', shell('make html'), delay=1)
    server.watch('_static/*', shell('make html'), delay=1)
    server.watch('_templates/*', shell('make html'), delay=1)

    server.serve(root='build/html', port=get_open_port())

def get_open_port():
    import socket
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("", 0))
    s.listen(1)
    port = s.getsockname()[1]
    s.close()
    return port

if __name__ == '__main__':
    main()