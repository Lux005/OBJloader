#ifndef  _LOADOBJ_H_  
#define  _LOADOBJ_H_  

#include <math.h> 
#include <iostream>  
#include <fstream>  
#include <string>  
#include <vector>  
#include "float3.h"
namespace OLoader{
	using namespace std;

	#define dot(u,v)   ((u).x * (v).x + (u).y * (v).y + (u).z*(v).z)    
	#define norm2(v)   dot(v,v)        // norm2 = squared length of vector    
	#define norm(v)    sqrt(norm2(v))  // norm = length of vector    
	#define d(u,v)     norm(u-v)       // distance = norm of difference    

	// Obj loader  
	struct TriangleFace
	{
		int v[3]; // vertex indices
		int n[3]; // vertex indices
		int t[3]; // vertex indices
	};


	struct TriangleMesh
	{
		vector<float3> verts;
		vector<TriangleFace> faces;
		vector<float3> vnormals;
		vector<float3> vtextures;
		bool textureON;
		string texturePath;
	  
		float3 bounding_box[2];
	
		float3 bounding_sphere_c;
		float bounding_sphere_r;
	};

	TriangleMesh mesh;

	int total_number_of_triangles = 0;

	// Scene bounding box   
	float3 scene_aabbox_min;
	float3 scene_aabbox_max;

	void loadObj(const std::string filename, TriangleMesh &mesh, bool MtlOn);
	void loadMtl(const std::string filename, TriangleMesh &mesh);

	void loadObj(const std::string filename, TriangleMesh &mesh, bool MtlOn=false)
	{
		std::ifstream in(filename.c_str());

		if (!in.good())
		{
			cout << "ERROR: loading obj:(" << filename << ") file is not good" << "\n";
			system("PAUSE");
			exit(0);
		}

		char buffer[256], str[255];
		float f1, f2, f3;
		mesh.textureON = false;
		while (!in.getline(buffer, 255).eof())
		{
			buffer[255] = '\0';

			sscanf_s(buffer, "%s", str, 255);
			
			if (MtlOn)
			{
				//reading mtllib
				if (memcmp(buffer, string("mtllib").c_str(), 6) == 0 && (buffer[6] == ' ' || buffer[6] == 32))
				{
					char mtlFile[255];
					if (sscanf_s(buffer, "mtllib %200[^\n]", &mtlFile, 255) == 1)
					{

						string mtlFname = mtlFile;
						loadMtl(mtlFname, mesh);
						continue;
					}
					cout << "ERROR: 'mtllib not in wanted format in OBJLoader'" << "\n";
					system("PAUSE");
					exit(-1);
				}
			}
			if (memcmp(buffer, string("texture").c_str(), 7) == 0 && (buffer[7] == ' ' || buffer[7] == 32))
			{
				char mtlFile[255];
				if (sscanf_s(buffer, "texture %200[^\n]", &mtlFile, 255) == 1)
				{
					mesh.textureON = true;
					mesh.texturePath = mtlFile;
					continue;
				}
				cout << "ERROR: 'texture' not in wanted format in OBJLoader'" << "\n";
				system("PAUSE");
				exit(-1);
			}
			// reading a vertex  
			if (buffer[0] == 'v' && (buffer[1] == ' ' || buffer[1] == 32))
			{
				if (sscanf(buffer, "v %f %f %f", &f1, &f2, &f3) == 3)
				{
					mesh.verts.push_back(float3(f1, f2, f3));
				}
				else
				{
					cout << "ERROR: 'v' not in wanted format in OBJLoader" << "\n";
					system("PAUSE");
					exit(-1);
				}
				continue;
			}
			if (buffer[0] == 'v' && buffer[1]=='t' && (buffer[2] == ' ' || buffer[2] == 32))
			{
				if (sscanf(buffer, "vt %f %f %f", &f1, &f2, &f3) == 3)
				{
					mesh.vtextures.push_back(float3(f1, f2, f3));
					continue;
				}
				if (sscanf(buffer, "vt %f %f", &f1, &f2) == 2)
				{
					mesh.vtextures.push_back(float3(f1, f2, 0));
					continue;
				}
				cout << "ERROR: 'vt' not in wanted format in OBJLoader" << "\n";
				system("PAUSE");
				exit(-1);
			}
			if (buffer[0] == 'v' && buffer[1] == 'n' && (buffer[2] == ' ' || buffer[2] == 32))
			{
				if (sscanf(buffer, "vn %f %f %f", &f1, &f2, &f3) == 3)
				{
					mesh.vnormals.push_back(float3(f1, f2, f3));
					continue;
				}
				cout << "ERROR: 'vn' not in wanted format in OBJLoader" << "\n";
				system("PAUSE");
				exit(-1);
			}
			// reading FaceMtls   
			if (buffer[0] == 'f' && (buffer[1] == ' ' || buffer[1] == 32))
			{
				TriangleFace f;
				int n4=0, t4=0,v4=0;
				for (int i = 0; i < 3; i++)
				{
					f.n[i] = 0;
					f.v[i] = 0;
					f.t[i] = 0;
				}
				int nt = sscanf(buffer, "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d", &f.v[0], &f.t[0], &f.n[0], &f.v[1], &f.t[1], &f.n[1], &f.v[2], &f.t[2], &f.n[2], &v4, &t4, &n4);
				if (nt == 12)
				{
					TriangleFace Tp;
					Tp.v[0] = f.v[2];
					Tp.v[1] = v4;
					Tp.v[2] = f.v[0];
					Tp.t[0] = f.t[2];
					Tp.t[1] = t4;
					Tp.t[2] = f.t[0];
					Tp.n[0] = f.n[2];
					Tp.n[1] = n4;
					Tp.n[2] = f.n[0];

					mesh.faces.push_back(Tp);
					mesh.faces.push_back(f);
					continue;
				}
				nt = sscanf(buffer, "f %d//%d %d//%d %d//%d %d//%d", &f.v[0], &f.n[0], &f.v[1], &f.n[1], &f.v[2], &f.n[2], &v4, &n4);
				if (nt == 8)
				{
					TriangleFace Tp;
					Tp.v[0] = f.v[2];
					Tp.v[1] = v4;
					Tp.v[2] = f.v[0];
					Tp.t[0] = f.t[2];
					Tp.t[1] = t4;
					Tp.t[2] = f.t[0];
					Tp.n[0] = f.n[2];
					Tp.n[1] = n4;
					Tp.n[2] = f.n[0];

					mesh.faces.push_back(Tp);
					mesh.faces.push_back(f);
					continue;
				}
				nt = sscanf(buffer, "f %d/%d %d/%d %d/%d %d/%d", &f.v[0], &f.t[0], &f.v[1], &f.t[1], &f.v[2], &f.t[2], &v4, &t4);
				if (nt == 8)
				{
					TriangleFace Tp;
					Tp.v[0] = f.v[2];
					Tp.v[1] = v4;
					Tp.v[2] = f.v[0];
					Tp.t[0] = f.t[2];
					Tp.t[1] = t4;
					Tp.t[2] = f.t[0];
					Tp.n[0] = f.n[2];
					Tp.n[1] = n4;
					Tp.n[2] = f.n[0];

					mesh.faces.push_back(Tp);
					mesh.faces.push_back(f);
					continue;
				}
				nt = sscanf(buffer, "f %d %d %d %d", &f.v[0],  &f.v[1],  &f.v[2], &v4);
				if (nt == 4)
				{
					TriangleFace Tp;
					Tp.v[0] = f.v[2];
					Tp.v[1] = v4;
					Tp.v[2] = f.v[0];
					Tp.t[0] = f.t[2];
					Tp.t[1] = t4;
					Tp.t[2] = f.t[0];
					Tp.n[0] = f.n[2];
					Tp.n[1] = n4;
					Tp.n[2] = f.n[0];

					mesh.faces.push_back(Tp);
					mesh.faces.push_back(f);
					continue;
				}
				nt = sscanf(buffer, "f %d/%d/%d %d/%d/%d %d/%d/%d", &f.v[0], &f.t[0], &f.n[0], &f.v[1], &f.t[1], &f.n[1], &f.v[2], &f.t[2], &f.n[2]);
				if (nt == 9)
				{
					mesh.faces.push_back(f);
					continue;
				}
				nt = sscanf(buffer, "f %d//%d %d//%d %d//%d", &f.v[0], &f.n[0], &f.v[1],  &f.n[1], &f.v[2],  &f.n[2]);
				if (nt == 6)
				{
					mesh.faces.push_back(f);
					continue;
				}
				nt = sscanf(buffer, "f %d/%d %d/%d %d/%d", &f.v[0], &f.t[0], &f.v[1], &f.t[1], &f.v[2], &f.t[2]);
				if (nt == 6)
				{
					mesh.faces.push_back(f);
					continue;
				}
				nt = sscanf(buffer, "f %d %d %d", &f.v[0], &f.v[1],  &f.v[2]);
				if (nt == 3)
				{
					mesh.faces.push_back(f);
					continue;
				}
				cout << "ERROR: 'f' not in wanted format in OBJLoader" << "\n";
				system("PAUSE");
				exit(-1);
			}
		}

		float xmin, ymin, zmin, xmax, ymax, zmax;
		int Pxmin, Pxmax, Pymin, Pymax, Pzmin, Pzmax;

		//calculate the bounding sphere  
		xmin = xmax = mesh.verts[0].x;
		ymin = ymax = mesh.verts[0].y;
		zmin = zmax = mesh.verts[0].z;
		Pxmin = Pxmax = Pymin = Pymax = Pzmin = Pzmax = 0;

		//calculate the bounding box  
		mesh.bounding_box[0] = float3(mesh.verts[0].x, mesh.verts[0].y, mesh.verts[0].z);
		mesh.bounding_box[1] = float3(mesh.verts[0].x, mesh.verts[0].y, mesh.verts[0].z);

		for (unsigned int i = 1; i < mesh.verts.size(); i++)
		{
			//update min value  
			mesh.bounding_box[0].x = fmin(mesh.verts[i].x, mesh.bounding_box[0].x);
			mesh.bounding_box[0].y = fmin(mesh.verts[i].y, mesh.bounding_box[0].y);
			mesh.bounding_box[0].z = fmin(mesh.verts[i].z, mesh.bounding_box[0].z);

			//update max value  
			mesh.bounding_box[1].x = fmax(mesh.verts[i].x, mesh.bounding_box[1].x);
			mesh.bounding_box[1].y = fmax(mesh.verts[i].y, mesh.bounding_box[1].y);
			mesh.bounding_box[1].z = fmax(mesh.verts[i].z, mesh.bounding_box[1].z);

			//update the  x min and max  
			if (mesh.verts[i].x < xmin){
				xmin = mesh.verts[i].x;
				Pxmin = i;
			}
			else if (mesh.verts[i].x > xmax){
				xmax = mesh.verts[i].x;
				Pxmax = i;
			}
			//update the y min and max  
			if (mesh.verts[i].y < ymin){
				ymin = mesh.verts[i].y;
				Pymin = i;
			}
			else if (mesh.verts[i].y > ymax){
				ymax = mesh.verts[i].y;
				Pymax = i;
			}
			//update the z min and max  
			if (mesh.verts[i].z < zmin){
				zmin = mesh.verts[i].z;
				Pzmin = i;
			}
			else if (mesh.verts[i].z > zmax){
				zmax = mesh.verts[i].z;
				Pzmax = i;
			}
		}

		//calculate the bounding sphere  
		float3 dVx = mesh.verts[Pxmax] - mesh.verts[Pxmin];
		float3 dVy = mesh.verts[Pymax] - mesh.verts[Pymin];
		float3 dVz = mesh.verts[Pzmax] - mesh.verts[Pzmin];
		float dx2 = norm2(dVx);
		float dy2 = norm2(dVy);
		float dz2 = norm2(dVz);

		float3 center;
		float  radius2;
		float  radius;

		if (dx2 >= dy2 && dx2 >= dz2) {                 // x direction is largest extent    
			center = mesh.verts[Pxmin] + (dVx / 2.0);   // Center = midpoint of extremes    
			radius2 = norm2(mesh.verts[Pxmax] - center);// radius squared    
		}
		else if (dy2 >= dx2  && dy2 >= dz2){                // y direction is largest extent    
			center = mesh.verts[Pymin] + (dVy / 2.0);   // Center = midpoint of extremes    
			radius2 = norm2(mesh.verts[Pymax] - center);// radius squared    
		}
		else{
			center = mesh.verts[Pzmin] + (dVz / 2.0);   // Center = midpoint of extremes    
			radius2 = norm2(mesh.verts[Pzmax] - center);// radius squared     
		}

		radius = sqrt(radius2);

		// now check that all points V[i] are in the ball    
		// and if not, expand the ball just enough to include them    
		float3 dV;
		float dist2, dist;
		for (unsigned int i = 0; i<mesh.verts.size(); i++)
		{
			dV = mesh.verts[i] - center;
			dist2 = norm2(dV);
			if (dist2 <= radius2) // V[i] is inside the ball already    
				continue;

			// V[i] not in ball, so expand ball to include it    
			dist = sqrt(dist2);
			radius = (radius + dist) / 2.0;         // enlarge radius just enough    
			radius2 = radius * radius;
			center = center + ((dist - radius) / dist) * dV;   // shift Center toward V[i]    

		}

		mesh.bounding_sphere_c = center;
		mesh.bounding_sphere_r = radius;

		cout << "----------obj file loaded-------------" << endl;
		cout << "number of faces:" << mesh.faces.size() << " number of vertices:" << mesh.verts.size() << endl;
		cout << "obj bounding box: min:("
			<< mesh.bounding_box[0].x << "," << mesh.bounding_box[0].y << "," << mesh.bounding_box[0].z << ") max:("
			<< mesh.bounding_box[1].x << "," << mesh.bounding_box[1].y << "," << mesh.bounding_box[1].z << ")" << endl
			<< "obj bounding sphere center:" << mesh.bounding_sphere_c.x << "," << mesh.bounding_sphere_c.y << "," << mesh.bounding_sphere_c.z << endl
			<< "obj bounding sphere radius:" << mesh.bounding_sphere_r << endl;

	}
	void loadMtl(const std::string filename, TriangleMesh &mesh)
	{
		std::ifstream in(filename.c_str());
		if (!in.good())
		{
			cout << "ERROR: loading mtl:(" << filename << ") file is not good" << "\n";
			system("PAUSE");
			exit(0);
		}

		char buffer[256], str[255];

		while (!in.getline(buffer, 255).eof())
		{
			buffer[255] = '\0';

			sscanf_s(buffer, "%s", str, 255);
			if (memcmp(buffer, string("newmtl").c_str(), 6) == 0 && (buffer[6] == ' ' || buffer[6] == 32))
			{
				char mtlName[255];
				if (sscanf_s(buffer, "newmtl %s", &mtlName, 255) == 1)
				{
					continue;
				}
				cout << "ERROR: 'newmtl' not in wanted format in OBJLoader'" << "\n";
				system("PAUSE");
				exit(-1);
			}
			if (memcmp(buffer, string("newmtl").c_str(), 6) == 0 && (buffer[6] == ' ' || buffer[6] == 32))
			{
				char mtlName[255];
				if (sscanf_s(buffer, "newmtl %s", &mtlName, 255) == 1)
				{
					continue;
				}
				cout << "ERROR: 'newmtl' not in wanted format in OBJLoader'" << "\n";
				system("PAUSE");
				exit(-1);
			}
			// reading a vertex  
		}
	}
}
#endif  
