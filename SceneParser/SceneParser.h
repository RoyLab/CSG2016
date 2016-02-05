#ifndef _SCENE_PARSER_H
#define _SCENE_PARSER_H
#include <vector>
#include <typedefs.h>

namespace GS
{
	class BaseMesh;
	class Light;
	class Camera;
}
class SceneParser{

public:
	SceneParser() {}
	virtual ~SceneParser() {}
	virtual bool Initialize() = 0;
	virtual bool LoadScene(const char* pFilename) = 0; 
	virtual void ProcessScene(std::vector<GS::BaseMesh*>& meshs, std::vector<GS::Light*>& lights,GS::Camera& camera) = 0;
	virtual void Release() = 0; 
	static SceneParser* GetParser(const char* pFilename);

};



#endif 