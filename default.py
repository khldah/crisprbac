#!C:\Python39\python.exe
import cgi, cgitb

form = cgi.FieldStorage()

print("Content-type:text/html \r\n\r\n")
print("Hello World")
print("""
<!doctype html>
<html>
    <head>
        <title>CRISPRBAC</title>
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />

    </head>
    <body>
        Hello World!
    </body>
</html>
""")