flex = /usr/local/Cellar/flex/2.6.2/bin/flex
bison = /usr/local/Cellar/bison/3.0.4/bin/bison

#SceneParser: SceneParser.tab.cc SceneLexer.yy.c
#	g++ -std=gnu++11 -o SceneParser SceneParser.tab.cc SceneLexer.yy.cpp main.cpp

SceneLexer.yy.c: SceneLexer.l SceneParser.tab.cc
	$(flex) -o SceneLexer.yy.cpp SceneLexer.l

SceneParser.tab.cc: SceneParser.yy
	$(bison) SceneParser.yy

clean:
	rm SceneParser.tab.* SceneLexer.yy.* location.hh position.hh stack.hh SceneParser
