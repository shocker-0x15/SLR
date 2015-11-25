#ifndef SLRSceneGraph_SceneParsingDriver_hpp
#define SLRSceneGraph_SceneParsingDriver_hpp

#include <libSLR/defines.h>
#include "references.h"
#include "SceneParser.tab.hh"
#include "Scene.h"

#define YY_DECL \
SLRSceneGraph::SceneParser::symbol_type yylex(SLRSceneGraph::SceneParsingDriver &driver)
YY_DECL;

namespace SLRSceneGraph {
	struct SceneParsingDriver {
		Scene scene;
        RenderingContext context;
        
		std::string file;
		bool traceScanning;
		bool traceParsing;
		location currentLocation;
		std::map<std::string, Element> variables;

		SceneParsingDriver() : 
		traceScanning(false), traceParsing(false) {

		};
        ~SceneParsingDriver() { };

		void beginScan();
		void endScan();

		int32_t parse(const std::string &f, Scene* outScene, RenderingContext* outContext) {
            variables["root"] = Element(Type::Node, scene.rootNode());
            
            file = f;
            beginScan();
            SceneParser parser(*this);
            parser.set_debug_level(traceParsing);
            int res = parser.parse();
            endScan();
            
            *outScene = scene;
            *outContext = std::move(context);
            
            return res;
		};

		void error(const location &loc, const std::string &msg) {
			std::cerr << loc << ": " << msg << std::endl;
		};
		void error(const std::string &msg) {
			std::cerr << msg << std::endl;
		};
	};
}

#endif
