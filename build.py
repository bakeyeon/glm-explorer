# engine.R을 app.src.html에 넣어 배포용 index.html 하나로 만든다
s = open("app.src.html", encoding="utf8").read()
e = open("engine.R", encoding="utf8").read() + "\n" + open("py.R", encoding="utf8").read()
assert "</script" not in e
open("index.html", "w", encoding="utf8", newline="\n").write(s.replace("/*ENGINE*/", e))
