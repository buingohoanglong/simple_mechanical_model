#include <iostream>
#include <math.h>
#include <GL/glut.h>

using namespace std;


const float PI = 3.1415926;

const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 600;


class Color {
public:
    float red;
    float green;
    float blue;
    Color(float red, float green, float blue) {
        this->red = red;
        this->green = green;
        this->blue = blue;
    }

    Color(const Color& color) {
        red = color.red;
        green = color.green;
        blue = color.blue;
    }
};

const int COLORNUM = 14;

Color ColorArr[COLORNUM] = {
    Color(1.0, 0.0, 0.0), 
    Color(0.0, 1.0, 0.0),
    Color(0.0,  0.0, 1.0),
    Color(1.0, 1.0,  0.0),
    Color(1.0, 0.0, 1.0),
    Color(0.0, 1.0, 1.0),
    Color(0.3, 0.3, 0.3),
    Color(0.5, 0.5, 0.5),
    Color(0.9,  0.9, 0.9),
    Color(1.0, 0.5,  0.5), 
    Color(0.5, 1.0, 0.5), 
    Color(0.5, 0.5, 1.0),
    Color(0.0, 0.0, 0.0), 
    Color(1.0, 1.0, 1.0)
};

class Point3f {
public:
    float x, y, z;

    void set(float dx, float dy, float dz) {
        x = dx;
        y = dy;
        z = dz;
    }

    void set(Point3f &p)
    {
        x = p.x;
        y = p.y;
        z = p.z;
    }

    Point3f() { x = y = z = 0; }

    Point3f(float dx, float dy, float dz) {
        x = dx;
        y = dy;
        z = dz;
    }
};

class Vector3f
{
public:
    float x, y, z;
    void set(float dx, float dy, float dz) {
        x = dx;
        y = dy;
        z = dz;
    }

    void set(Vector3f &v) {
        x = v.x;
        y = v.y;
        z = v.z;
    }

    void flip() {
        x = -x;
        y = -y;
        z = -z;
    }

    void normalize() {
        float cardinality = sqrt(x*x + y*y + z*z);
        x = x / cardinality;
        y = y / cardinality;
        z = z / cardinality;
    }

    Vector3f() { x = y = z = 0; }

    Vector3f(float dx, float dy, float dz) {
        x = dx;
        y = dy;
        z = dz;
    }

    Vector3f(Vector3f &v) {
        x = v.x;
        y = v.y;
        z = v.z;
    }

    Vector3f cross(Vector3f b) {
        Vector3f v(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x);
        return v;
    }

    float dot(Vector3f b) {
        return x*b.x + y*b.y + z*b.z;
    }
};


class VertexID {
public:
	int	vertIndex;
	int	colorIndex;
};

class Face {
public:
	int	nVerts;
	VertexID* vert;
    Vector3f facenorm;

	Face() {
		nVerts	= 0;
		vert	= NULL;
        facenorm = Vector3f();
	}

	~Face() {
		if(vert !=NULL)
		{
			delete[] vert;
			vert = NULL;
		}
		nVerts = 0;
	}
};

class Mesh {
public:
	int numVerts;
	Point3f* pt;
	int	numFaces;
	Face* face;
    Color color = Color(0,0,0);
	GLfloat ambient[4] = {0.2, 0.2, 0.2, 1.0};
	GLfloat diffuse[4] = {0.8, 0.8, 0.8, 0.0};
	GLfloat specular[4] = {1.0, 1.0, 1.0, 0.0};
	GLfloat shiness = 50.0f;

public:
	Mesh() {
		numVerts	= 0;
		pt		= NULL;
		numFaces	= 0;
		face		= NULL;
	}

	~Mesh() {
		if (pt != NULL)
		{
			delete[] pt;
		}
		if(face != NULL)
		{
			delete[] face;
		}
		numVerts = 0;
		numFaces = 0;
	}

	void DrawWireframe() {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glColor3f(color.red, color.green, color.blue);
        for (int f = 0; f < numFaces; f++) {
            glBegin(GL_POLYGON);
            for (int v = 0; v < face[f].nVerts; v++)
            {
                int	iv = face[f].vert[v].vertIndex;
                glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
            }
            glEnd();
        }
    }
	
    void DrawColor() {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f(color.red, color.green, color.blue);
        for (int f = 0; f < numFaces; f++)
        {
            glBegin(GL_POLYGON);
            for (int v = 0; v < face[f].nVerts; v++)
            {
                int	iv = face[f].vert[v].vertIndex;
                // int	ic = face[f].vert[v].colorIndex; 
                glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
            }
            glEnd();
        }
    }

    void setColor(Color color) {
        this->color = color;
    }

    void Draw() {
        // set material
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, this->ambient);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, this->diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, this->specular);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, this->shiness);

        // draw
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        for (int f = 0; f < numFaces; f++){
            glBegin(GL_POLYGON);
            for (int v = 0; v < face[f].nVerts; v++){
                int	iv = face[f].vert[v].vertIndex;
                glNormal3f(face[f].facenorm.x, face[f].facenorm.y, face[f].facenorm.z); // use face norm
                glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
            }
            glEnd();
        }
    }

    void setupMaterial(Color color) {
        float newcolor[4];
        newcolor[0] = color.red;
        newcolor[1] = color.green;
        newcolor[2] = color.blue;
        newcolor[3] = 1.0;

        if(ambient != NULL) {
            for(int i=0; i<4; i++) this->ambient[i] = newcolor[i];
        }
        if(diffuse != NULL) {
            for(int i=0; i<4; i++) this->diffuse[i] = newcolor[i];
        }
        // if(diffuse != NULL) {
        //     for(int i=0; i<4; i++) this->specular[i] = specular[i];
        // }
        // this->shiness = shiness;
    }

    // only call after init vertices and faces
    // calculate face norms using Newell's method
    void CalculateFacesNorm() {
        for(int f=0; f<numFaces; f++) {
            float mx = 0;
            float my = 0;
            float mz = 0;
            for(int i=0; i<face[f].nVerts; i++) {
                int vertexIndex = face[f].vert[i].vertIndex;
                int vertexIndexNext = face[f].vert[(i+1) % face[f].nVerts].vertIndex;
                mx += (pt[vertexIndex].y - pt[vertexIndexNext].y) * (pt[vertexIndex].z + pt[vertexIndexNext].z);
                my += (pt[vertexIndex].z - pt[vertexIndexNext].z) * (pt[vertexIndex].x + pt[vertexIndexNext].x);
                mz += (pt[vertexIndex].x - pt[vertexIndexNext].x) * (pt[vertexIndex].y + pt[vertexIndexNext].y);
            }
            face[f].facenorm.set(mx, my, mz);
            face[f].facenorm.normalize();
        }
    }

    virtual void init() = 0;
};

class Tile: public Mesh {
public:
    Color firstColor = Color(0,0,0);
    Color secondColor = Color(0,0,0);
	// GLfloat ambient1[4] = {0.2, 0.2, 0.2, 1.0};
	// GLfloat diffuse1[4] = {0.8, 0.8, 0.8, 1.0};
	// GLfloat specular1[4] = {1.0, 1.0, 1.0, 1.0};
	// GLfloat ambient2[4] = {0.2, 0.2, 0.2, 1.0};
	// GLfloat diffuse2[4] = {0.8, 0.8, 0.8, 1.0};
	// GLfloat specular2[4] = {1.0, 1.0, 1.0, 1.0};
	// GLfloat shiness = 50.0f;
    Tile(): Mesh() {}
    ~Tile() {}

    void setupFirstColor(Color color) {
       firstColor = color;
    }

    void setupSecondColor(Color color) {
        secondColor = color;
    }

    void DrawColor() {
        glDisable(GL_LIGHTING);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f(secondColor.red, secondColor.green, secondColor.blue);
        for (int f = 0; f < numFaces - 1; f++)
        {
            glBegin(GL_POLYGON);
            for (int v = 0; v < face[f].nVerts; v++)
            {
                int	iv = face[f].vert[v].vertIndex;
                // int	ic = face[f].vert[v].colorIndex; 
                glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
            }
            glEnd();
        }

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f(firstColor.red, firstColor.green, firstColor.blue);
        glBegin(GL_POLYGON);
        for (int v = 0; v < face[numFaces-1].nVerts; v++)
        {
            int	iv = face[numFaces-1].vert[v].vertIndex;
            // int	ic = face[f].vert[v].colorIndex; 
            glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
        }
        glEnd();

        glEnable(GL_LIGHTING);
    }
};

class RectangularPrism : public Mesh {
public:
    float width, length, height;

    RectangularPrism(float width, float length, float height): Mesh(), width(width), length(length), height(height) {}
    ~RectangularPrism() {}

    virtual void init() {
        // init vertices
        numVerts = 8;
        pt = new Point3f[numVerts];
        // bottom vertices
        pt[0].set(0, 0, 0);
        pt[1].set(0, 0, width);
        pt[2].set(length, 0, width);
        pt[3].set(length, 0, 0);
        // top vertices
        pt[4].set(0, height, 0);
        pt[5].set(0, height, width);
        pt[6].set(length, height, width);
        pt[7].set(length, height, 0);

        // init faces
        numFaces = 6;
        face = new Face[numFaces];
        // bottom
        face[0].nVerts = 4;
        face[0].vert = new VertexID[face[0].nVerts];
        face[0].vert[0].vertIndex = 0;
        face[0].vert[1].vertIndex = 1;
        face[0].vert[2].vertIndex = 2;
        face[0].vert[3].vertIndex = 3;
        // top
        face[1].nVerts = 4;
        face[1].vert = new VertexID[face[1].nVerts];
        face[1].vert[0].vertIndex = 4;
        face[1].vert[1].vertIndex = 7;
        face[1].vert[2].vertIndex = 6;
        face[1].vert[3].vertIndex = 5;
        // left
        face[2].nVerts = 4;
        face[2].vert = new VertexID[face[2].nVerts];
        face[2].vert[0].vertIndex = 4;
        face[2].vert[1].vertIndex = 5;
        face[2].vert[2].vertIndex = 1;
        face[2].vert[3].vertIndex = 0;
        // right
        face[3].nVerts = 4;
        face[3].vert = new VertexID[face[3].nVerts];
        face[3].vert[0].vertIndex = 6;
        face[3].vert[1].vertIndex = 7;
        face[3].vert[2].vertIndex = 3;
        face[3].vert[3].vertIndex = 2;
        // front
        face[4].nVerts = 4;
        face[4].vert = new VertexID[face[4].nVerts];
        face[4].vert[0].vertIndex = 1;
        face[4].vert[1].vertIndex = 5;
        face[4].vert[2].vertIndex = 6;
        face[4].vert[3].vertIndex = 2;
        // back
        face[5].nVerts = 4;
        face[5].vert = new VertexID[face[5].nVerts];
        face[5].vert[0].vertIndex = 4;
        face[5].vert[1].vertIndex = 0;
        face[5].vert[2].vertIndex = 3;
        face[5].vert[3].vertIndex = 7;
            
        CalculateFacesNorm();
    }
};

class Cylinder : public Mesh {
public:
    float radius, height;
    int nSegments;

    Cylinder(float radius, float height, int nSegments): Mesh(), radius(radius), height(height), nSegments(nSegments) {}
    ~Cylinder() {}

    virtual void init() {
        // init vertices
        numVerts = nSegments * 2 + 2;
        pt = new Point3f[numVerts];
        // two center
        pt[0].set(0, 0, 0);
        pt[1].set(0, height, 0);

        // vertices in 2 circles
        float angle = 0;
        float deltaAngle = 2*PI / nSegments;
        for(int i=0; i<nSegments; i++) {
            // bottom vertex
            pt[i + 2].set(radius*sin(angle), 0, radius*cos(angle));
            // top vertex
            pt[i + 2 + nSegments].set(radius*sin(angle), height, radius*cos(angle));
            angle += deltaAngle;
        }

        // init faces
        numFaces = 3*nSegments;
        face = new Face[numFaces];
        for(int i=0; i<nSegments; i++) {
            // bottom faces
            face[i].nVerts = 3;
            face[i].vert = new VertexID[face[i].nVerts];
            face[i].vert[0].vertIndex = 0;
            face[i].vert[1].vertIndex = i % nSegments + 2;
            face[i].vert[2].vertIndex = (i+1) % nSegments + 2;

            // top faces
            face[i + nSegments].nVerts = 3;
            face[i + nSegments].vert = new VertexID[face[i + nSegments].nVerts];
            face[i + nSegments].vert[0].vertIndex = 1;
            face[i + nSegments].vert[1].vertIndex = (i+1) % nSegments + 2 + nSegments;
            face[i + nSegments].vert[2].vertIndex = i % nSegments + 2 + nSegments;
            

            // side faces
            face[i + 2*nSegments].nVerts = 4;
            face[i + 2*nSegments].vert = new VertexID[face[i + 2*nSegments].nVerts];
            face[i + 2*nSegments].vert[0].vertIndex = i % nSegments + 2;
            face[i + 2*nSegments].vert[1].vertIndex = i % nSegments + 2 + nSegments;
            face[i + 2*nSegments].vert[2].vertIndex = (i+1) % nSegments + 2 + nSegments;   
            face[i + 2*nSegments].vert[3].vertIndex = (i+1) % nSegments + 2;     
        }

        CalculateFacesNorm();
            
    }
};

class Crank: public Mesh {
public:
    float width, height, length;
    int nSegments;
    Crank(float width, float height, float length, int nSegments): Mesh(),
        width(width), height(height), length(length), nSegments(nSegments) {}
    ~Crank() {}

    virtual void init() {
        // init vertices;
        numVerts = (4 + 2 + 2*(nSegments+1)) *2;
        pt = new Point3f[numVerts];

        // rectangular prism vertices
        // bottom vertices
        pt[numVerts - 1].set(-(length - width) / 2, 0, -width / 2);
        pt[numVerts - 2].set((length - width) / 2, 0, -width / 2);
        pt[numVerts - 3].set((length - width) / 2, 0, width / 2);
        pt[numVerts - 4].set(-(length - width) / 2, 0, width / 2);
        // top vertices
        pt[numVerts - 5].set(-(length - width) / 2, height, -width / 2);
        pt[numVerts - 6].set((length - width) / 2, height, -width / 2);
        pt[numVerts - 7].set((length - width) / 2, height, width / 2);
        pt[numVerts - 8].set(-(length - width) / 2, height, width / 2);

        // centers
        pt[numVerts - 9].set(-(length - width) / 2, 0, 0);
        pt[numVerts - 10].set(-(length - width) / 2, height, 0);
        pt[numVerts - 11].set((length - width) / 2, 0, 0);
        pt[numVerts - 12].set((length - width) / 2, height, 0);

        // circle vertices
        float radius = width / 2;
        float angle = 0;
        float deltaAngle = PI / nSegments;
        for(int i=0; i<nSegments + 1; i++) {
            // bottom left
            pt[i].set(-(length - width) / 2 - radius*sin(angle), 0, radius*cos(angle));
            // top left
            pt[i + nSegments + 1].set(-(length - width) / 2 - radius*sin(angle), height, radius*cos(angle));
            // bottom right
            pt[i + 2*(nSegments + 1)].set((length - width) / 2 + radius*sin(angle), 0, radius*cos(angle));
            // top right
            pt[i + 3*(nSegments + 1)].set((length - width) / 2 + radius*sin(angle), height, radius*cos(angle));
            angle += deltaAngle;
        }


        // init faces
        numFaces = 4 + (3*nSegments)*2;
        face = new Face[numFaces];

        // rectangular faces
        // bottom
        face[numFaces - 1].nVerts = 4;
        face[numFaces - 1].vert = new VertexID[face[numFaces - 1].nVerts];
        face[numFaces - 1].vert[0].vertIndex = numVerts - 1;
        face[numFaces - 1].vert[1].vertIndex = numVerts - 4;
        face[numFaces - 1].vert[2].vertIndex = numVerts - 3;
        face[numFaces - 1].vert[3].vertIndex = numVerts - 2;
        // top
        face[numFaces - 2].nVerts = 4;
        face[numFaces - 2].vert = new VertexID[face[numFaces - 2].nVerts];
        face[numFaces - 2].vert[0].vertIndex = numVerts - 5;
        face[numFaces - 2].vert[1].vertIndex = numVerts - 6;
        face[numFaces - 2].vert[2].vertIndex = numVerts - 7;
        face[numFaces - 2].vert[3].vertIndex = numVerts - 8;
        // front
        face[numFaces - 3].nVerts = 4;
        face[numFaces - 3].vert = new VertexID[face[numFaces - 3].nVerts];
        face[numFaces - 3].vert[0].vertIndex = numVerts - 4;
        face[numFaces - 3].vert[1].vertIndex = numVerts - 8;
        face[numFaces - 3].vert[2].vertIndex = numVerts - 7;
        face[numFaces - 3].vert[3].vertIndex = numVerts - 3;
        // back
        face[numFaces - 4].nVerts = 4;
        face[numFaces - 4].vert = new VertexID[face[numFaces - 4].nVerts];
        face[numFaces - 4].vert[0].vertIndex = numVerts - 1;
        face[numFaces - 4].vert[1].vertIndex = numVerts - 2;
        face[numFaces - 4].vert[2].vertIndex = numVerts - 6;
        face[numFaces - 4].vert[3].vertIndex = numVerts - 5;

        // cylinder faces
        for(int i=0; i<nSegments; i++) {
            // bottom left
            face[i].nVerts = 3;
            face[i].vert = new VertexID[face[i].nVerts];
            face[i].vert[0].vertIndex = numVerts - 9;
            face[i].vert[1].vertIndex = i+1;
            face[i].vert[2].vertIndex = i;
            // top left
            face[i + nSegments].nVerts = 3;
            face[i + nSegments].vert = new VertexID[face[i + nSegments].nVerts];
            face[i + nSegments].vert[0].vertIndex = numVerts - 10;
            face[i + nSegments].vert[1].vertIndex = i + nSegments+1;
            face[i + nSegments].vert[2].vertIndex = (i+1) + nSegments+1;
            // bottom right
            face[i + 2*nSegments].nVerts = 3;
            face[i + 2*nSegments].vert = new VertexID[face[i + 2*nSegments].nVerts];
            face[i + 2*nSegments].vert[0].vertIndex = numVerts - 11;
            face[i + 2*nSegments].vert[1].vertIndex = i  + 2*(nSegments+1);
            face[i + 2*nSegments].vert[2].vertIndex = (i+1)  + 2*(nSegments+1);
            // top right
            face[i + 3*nSegments].nVerts = 3;
            face[i + 3*nSegments].vert = new VertexID[face[i + 3*nSegments].nVerts];
            face[i + 3*nSegments].vert[0].vertIndex = numVerts - 12;
            face[i + 3*nSegments].vert[1].vertIndex = (i+1)  + 3*(nSegments+1);
            face[i + 3*nSegments].vert[2].vertIndex = i  + 3*(nSegments+1);

            // side left
            face[i + 4*nSegments].nVerts = 4;
            face[i + 4*nSegments].vert = new VertexID[face[i + 4*nSegments].nVerts];
            face[i + 4*nSegments].vert[0].vertIndex = (i+1) ;
            face[i + 4*nSegments].vert[1].vertIndex = (i+1)  + nSegments+1;
            face[i + 4*nSegments].vert[2].vertIndex = i  + nSegments+1;    
            face[i + 4*nSegments].vert[3].vertIndex = i ;        
            // side right
            face[i + 5*nSegments].nVerts = 4;
            face[i + 5*nSegments].vert = new VertexID[face[i + 5*nSegments].nVerts];
            face[i + 5*nSegments].vert[0].vertIndex = i  + 2*(nSegments+1);
            face[i + 5*nSegments].vert[1].vertIndex = i  + 3*(nSegments+1);
            face[i + 5*nSegments].vert[2].vertIndex = (i+1)  + 3*(nSegments+1);   
            face[i + 5*nSegments].vert[3].vertIndex = (i+1)  + 2*(nSegments+1);  
        }

        CalculateFacesNorm();

    }
};

class Bearing : public Mesh {
public:
    float bottomWidth, bottomHeight, bottomLength, bodyLength, bodyHeight, radius;
    int nSegments;
    Bearing(float bottomWidth, float bottomHeight, float bottomLength, 
        float bodyLength, float bodyHeight, float radius, int nSegments):
            Mesh(), bottomWidth(bottomWidth), bottomHeight(bottomHeight), bottomLength(bottomLength),
            bodyLength(bodyLength), bodyHeight(bodyHeight), radius(radius), nSegments(nSegments) {}
    ~Bearing() {}

    virtual void init() {
        // init vertices
        numVerts = (12 + 2*nSegments)*2;
        pt = new Point3f[numVerts];

        // rectangular bottom
        pt[numVerts - 1].set(-bottomLength/2, bottomHeight, bottomWidth/2);
        pt[numVerts - 2].set(-bottomLength/2, 0, bottomWidth/2);
        pt[numVerts - 3].set(bottomLength/2, 0, bottomWidth/2);
        pt[numVerts - 4].set(bottomLength/2, bottomHeight, bottomWidth/2);

        pt[numVerts - 13].set(-bottomLength/2, bottomHeight, -bottomWidth/2);
        pt[numVerts - 14].set(-bottomLength/2, 0, -bottomWidth/2);
        pt[numVerts - 15].set(bottomLength/2, 0, -bottomWidth/2);
        pt[numVerts - 16].set(bottomLength/2, bottomHeight, -bottomWidth/2);

        // rectangular body
        // left
        pt[numVerts - 5].set(-bodyLength/2, bottomHeight, bottomWidth/2);
        pt[numVerts - 6].set(-radius, bottomHeight, bottomWidth/2);
        pt[numVerts - 9].set(-bodyLength/2, bottomHeight + (bodyHeight - (bodyLength/2 - radius) - radius), bottomWidth/2);
        pt[numVerts - 10].set(-radius, bottomHeight + (bodyHeight - (bodyLength/2 - radius) - radius), bottomWidth/2);

        pt[numVerts - 17].set(-bodyLength/2, bottomHeight, -bottomWidth/2);
        pt[numVerts - 18].set(-radius, bottomHeight, -bottomWidth/2);
        pt[numVerts - 21].set(-bodyLength/2, bottomHeight + (bodyHeight - (bodyLength/2 - radius) - radius), -bottomWidth/2);
        pt[numVerts - 22].set(-radius, bottomHeight + (bodyHeight - (bodyLength/2 - radius) - radius), -bottomWidth/2);
        // right
        pt[numVerts - 7].set(radius, bottomHeight, bottomWidth/2);
        pt[numVerts - 8].set(bodyLength/2, bottomHeight, bottomWidth/2);
        pt[numVerts - 11].set(radius, bottomHeight + (bodyHeight - (bodyLength/2 - radius) - radius), bottomWidth/2);
        pt[numVerts - 12].set(bodyLength/2, bottomHeight + (bodyHeight - (bodyLength/2 - radius) - radius), bottomWidth/2);

        pt[numVerts - 19].set(radius, bottomHeight, -bottomWidth/2);
        pt[numVerts - 20].set(bodyLength/2, bottomHeight, -bottomWidth/2);
        pt[numVerts - 23].set(radius, bottomHeight + (bodyHeight - (bodyLength/2 - radius) - radius), -bottomWidth/2);
        pt[numVerts - 24].set(bodyLength/2, bottomHeight + (bodyHeight - (bodyLength/2 - radius) - radius), -bottomWidth/2);


        
        float angle = 0;
        float deltaAngle = 2*PI / nSegments;
        for(int i=0; i<nSegments; i++) {
            //inner circular vertices
            // front
            pt[i].set(radius*cos(angle), bottomHeight + bodyHeight - bodyLength/2 + radius*sin(angle), bottomWidth / 2);
            // back
            pt[i + nSegments].set(radius*cos(angle), bottomHeight + bodyHeight - bodyLength/2 + radius*sin(angle), -bottomWidth / 2);
            
            if(i < nSegments / 2) {
                // outer semi-circular vertices
                // front
                pt[i + 2*nSegments].set(bodyLength/2*cos(angle), bottomHeight + bodyHeight - bodyLength/2 + bodyLength/2*sin(angle), bottomWidth / 2);
                // back
                pt[i + 2*nSegments + nSegments / 2].set(bodyLength/2*cos(angle), bottomHeight + bodyHeight - bodyLength/2 + bodyLength/2*sin(angle), -bottomWidth / 2);
            } else {
                // inner circle - bottom connect face vertices
                // front
                pt[i + 3*nSegments - nSegments  /2].set(radius*cos(angle), bottomHeight, bottomWidth / 2);
                // back
                pt[i + 3*nSegments].set(radius*cos(angle), bottomHeight, -bottomWidth / 2);
            }

            angle += deltaAngle;
        }






        // init faces
        numFaces = 6 + 3*2 + nSegments*4;
        face = new Face[numFaces];

        // rectangular bottom
        // bottom
        face[numFaces - 1].nVerts = 4;
        face[numFaces - 1].vert = new VertexID[face[numFaces - 1].nVerts];
        face[numFaces - 1].vert[0].vertIndex = numVerts-2;
        face[numFaces - 1].vert[1].vertIndex = numVerts-3;
        face[numFaces - 1].vert[2].vertIndex = numVerts-15;
        face[numFaces - 1].vert[3].vertIndex = numVerts-14;
        // top
        face[numFaces - 2].nVerts = 4;
        face[numFaces - 2].vert = new VertexID[face[numFaces - 2].nVerts];
        face[numFaces - 2].vert[0].vertIndex = numVerts-1;
        face[numFaces - 2].vert[1].vertIndex = numVerts-13;
        face[numFaces - 2].vert[2].vertIndex = numVerts-16;
        face[numFaces - 2].vert[3].vertIndex = numVerts-4;
        // left
        face[numFaces - 3].nVerts = 4;
        face[numFaces - 3].vert = new VertexID[face[numFaces - 3].nVerts];
        face[numFaces - 3].vert[0].vertIndex = numVerts-1;
        face[numFaces - 3].vert[1].vertIndex = numVerts-2;
        face[numFaces - 3].vert[2].vertIndex = numVerts-14;
        face[numFaces - 3].vert[3].vertIndex = numVerts-13;
        // right
        face[numFaces - 4].nVerts = 4;
        face[numFaces - 4].vert = new VertexID[face[numFaces - 4].nVerts];
        face[numFaces - 4].vert[0].vertIndex = numVerts-3;
        face[numFaces - 4].vert[1].vertIndex = numVerts-4;
        face[numFaces - 4].vert[2].vertIndex = numVerts-16;
        face[numFaces - 4].vert[3].vertIndex = numVerts-15;
        // front
        face[numFaces - 5].nVerts = 4;
        face[numFaces - 5].vert = new VertexID[face[numFaces - 5].nVerts];
        face[numFaces - 5].vert[0].vertIndex = numVerts-1;
        face[numFaces - 5].vert[1].vertIndex = numVerts-4;
        face[numFaces - 5].vert[2].vertIndex = numVerts-3;
        face[numFaces - 5].vert[3].vertIndex = numVerts-2;
        // back
        face[numFaces - 6].nVerts = 4;
        face[numFaces - 6].vert = new VertexID[face[numFaces - 6].nVerts];
        face[numFaces - 6].vert[0].vertIndex = numVerts-13;
        face[numFaces - 6].vert[1].vertIndex = numVerts-14;
        face[numFaces - 6].vert[2].vertIndex = numVerts-15;
        face[numFaces - 6].vert[3].vertIndex = numVerts-16;


        // rectangular body left
        // front
        face[numFaces - 7].nVerts = 4;
        face[numFaces - 7].vert = new VertexID[face[numFaces - 7].nVerts];
        face[numFaces - 7].vert[0].vertIndex = numVerts-5;
        face[numFaces - 7].vert[1].vertIndex = numVerts-9;
        face[numFaces - 7].vert[2].vertIndex = numVerts-10;
        face[numFaces - 7].vert[3].vertIndex = numVerts-6;
        // left
        face[numFaces - 8].nVerts = 4;
        face[numFaces - 8].vert = new VertexID[face[numFaces - 8].nVerts];
        face[numFaces - 8].vert[0].vertIndex = numVerts-9;
        face[numFaces - 8].vert[1].vertIndex = numVerts-5;
        face[numFaces - 8].vert[2].vertIndex = numVerts-17;
        face[numFaces - 8].vert[3].vertIndex = numVerts-21;
        // back
        face[numFaces - 9].nVerts = 4;
        face[numFaces - 9].vert = new VertexID[face[numFaces - 9].nVerts];
        face[numFaces - 9].vert[0].vertIndex = numVerts-22;
        face[numFaces - 9].vert[1].vertIndex = numVerts-21;
        face[numFaces - 9].vert[2].vertIndex = numVerts-17;
        face[numFaces - 9].vert[3].vertIndex = numVerts-18;   

        // rectangular body right
        // front
        face[numFaces - 10].nVerts = 4;
        face[numFaces - 10].vert = new VertexID[face[numFaces - 10].nVerts];
        face[numFaces - 10].vert[0].vertIndex = numVerts-7;
        face[numFaces - 10].vert[1].vertIndex = numVerts-11;
        face[numFaces - 10].vert[2].vertIndex = numVerts-12;
        face[numFaces - 10].vert[3].vertIndex = numVerts-8;
        // right
        face[numFaces - 11].nVerts = 4;
        face[numFaces - 11].vert = new VertexID[face[numFaces - 11].nVerts];
        face[numFaces - 11].vert[0].vertIndex = numVerts-8;
        face[numFaces - 11].vert[1].vertIndex = numVerts-12;
        face[numFaces - 11].vert[2].vertIndex = numVerts-24;
        face[numFaces - 11].vert[3].vertIndex = numVerts-20;
        // back
        face[numFaces - 12].nVerts = 4;
        face[numFaces - 12].vert = new VertexID[face[numFaces - 12].nVerts];
        face[numFaces - 12].vert[0].vertIndex = numVerts-24;
        face[numFaces - 12].vert[1].vertIndex = numVerts-23;
        face[numFaces - 12].vert[2].vertIndex = numVerts-19;
        face[numFaces - 12].vert[3].vertIndex = numVerts-20;

        
        for(int i=0; i<nSegments; i++) {
            // inner circular faces
            face[i].nVerts = 4;
            face[i].vert = new VertexID[face[i].nVerts];
            face[i].vert[0].vertIndex = i % nSegments;
            face[i].vert[1].vertIndex = (i+1) % nSegments;
            face[i].vert[2].vertIndex = (i+1) % nSegments + nSegments;
            face[i].vert[3].vertIndex = i % nSegments + nSegments;

            if(i < nSegments / 2 ) {
                // inner-outer circular connect faces
                // front
                face[i + nSegments].nVerts = 4;
                face[i + nSegments].vert = new VertexID[face[i + nSegments].nVerts];
                bool isLast = (i == nSegments / 2 - 1);
                face[i + nSegments].vert[0].vertIndex =  i % nSegments;
                face[i + nSegments].vert[1].vertIndex = isLast ? numVerts - 10 : (i+1) % nSegments;
                face[i + nSegments].vert[2].vertIndex = isLast ? numVerts - 9: (i+1) % nSegments + 2*nSegments;
                face[i + nSegments].vert[3].vertIndex = i % nSegments + 2*nSegments; 
                // back
                face[i + nSegments + nSegments/2].nVerts = 4;
                face[i + nSegments + nSegments/2].vert = new VertexID[face[i + nSegments + nSegments/2].nVerts];
                face[i + nSegments + nSegments/2].vert[0].vertIndex = isLast ? numVerts - 22 : (i+1) % nSegments + nSegments;
                face[i + nSegments + nSegments/2].vert[1].vertIndex = i % nSegments + nSegments;
                face[i + nSegments + nSegments/2].vert[2].vertIndex = i % nSegments + 2*nSegments + nSegments/2;
                face[i + nSegments + nSegments/2].vert[3].vertIndex = isLast ? numVerts - 21 : (i+1) % nSegments + 2*nSegments + nSegments/2; 
            
                // outer semi-circular faces
                face[i + 3*nSegments].nVerts = 4;
                face[i + 3*nSegments].vert = new VertexID[face[i + 3*nSegments].nVerts];
                // bool isLast = (i == nSegments / 2 - 1);
                face[i + 3*nSegments].vert[0].vertIndex =  i % nSegments + 2*nSegments;
                face[i + 3*nSegments].vert[1].vertIndex = isLast ? numVerts - 9 : (i+1) % nSegments + 2*nSegments;
                face[i + 3*nSegments].vert[2].vertIndex = isLast ? numVerts - 21: (i+1) % nSegments + 2*nSegments + nSegments/2;
                face[i + 3*nSegments].vert[3].vertIndex = i % nSegments + 2*nSegments + nSegments/2; 
            } else {
                // inner circle - bottom connect faces
                bool isLast = (i == nSegments - 1);
                // front
                face[i + 2*nSegments - nSegments/2].nVerts = 4;
                face[i + 2*nSegments - nSegments/2].vert = new VertexID[face[i + 2*nSegments - nSegments/2].nVerts];
                face[i + 2*nSegments - nSegments/2].vert[0].vertIndex = i % nSegments + 3*nSegments - nSegments  /2;
                face[i + 2*nSegments - nSegments/2].vert[1].vertIndex = i % nSegments;
                face[i + 2*nSegments - nSegments/2].vert[2].vertIndex = isLast ? numVerts - 11 : (i+1) % nSegments;
                face[i + 2*nSegments - nSegments/2].vert[3].vertIndex = isLast ? numVerts - 7 : (i+1) % nSegments + 3*nSegments - nSegments  /2;
                // back
                face[i + 2*nSegments].nVerts = 4;
                face[i + 2*nSegments].vert = new VertexID[face[i + 2*nSegments].nVerts];
                face[i + 2*nSegments].vert[0].vertIndex = i % nSegments + nSegments;
                face[i + 2*nSegments].vert[1].vertIndex = i % nSegments + 3*nSegments;
                face[i + 2*nSegments].vert[2].vertIndex = isLast ? numVerts - 19 : (i+1) % nSegments + 3*nSegments;
                face[i + 2*nSegments].vert[3].vertIndex = isLast ? numVerts - 23 : (i+1) % nSegments + nSegments;
            }
           
        }


        CalculateFacesNorm();
    }
};

class UBar: public Mesh {
public:
    float width, height, length;
    UBar(float width, float height, float length): Mesh(), width(width), height(height), length(length) {}
    ~UBar() {}

    virtual void init() {
        // init vertices
        numVerts = 16;
        pt = new Point3f[numVerts];
        // front
        pt[0].set(-length/2, 0, width);
        pt[1].set(-length/2, height, width);
        pt[2].set(-length/2 + width, height, width);
        pt[3].set(-length/2 + width, width, width);
        pt[4].set(length/2 - width, width, width);
        pt[5].set(length/2 - width, height, width);
        pt[6].set(length/2, height, width);
        pt[7].set(length/2, 0, width);

        // back
        pt[8].set(-length/2, 0, 0);
        pt[9].set(-length/2, height, 0);
        pt[10].set(-length/2 + width, height, 0);
        pt[11].set(-length/2 + width, width, 0);
        pt[12].set(length/2 - width, width, 0);
        pt[13].set(length/2 - width, height, 0);
        pt[14].set(length/2, height, 0);
        pt[15].set(length/2, 0, 0);


        // init faces
        numFaces = 14;
        face = new Face[numFaces];

        // bottom bar
        // bottom
        face[0].nVerts = 4;
        face[0].vert = new VertexID[face[0].nVerts];
        face[0].vert[0].vertIndex = 0;
        face[0].vert[1].vertIndex = 7;
        face[0].vert[2].vertIndex = 15;
        face[0].vert[3].vertIndex = 8;
        // top
        face[1].nVerts = 4;
        face[1].vert = new VertexID[face[1].nVerts];
        face[1].vert[0].vertIndex = 3;
        face[1].vert[1].vertIndex = 11;
        face[1].vert[2].vertIndex = 12;
        face[1].vert[3].vertIndex = 4;
        // front
        face[2].nVerts = 4;
        face[2].vert = new VertexID[face[2].nVerts];
        face[2].vert[0].vertIndex = 0;
        face[2].vert[1].vertIndex = 3;
        face[2].vert[2].vertIndex = 4;
        face[2].vert[3].vertIndex = 7;
        // back
        face[3].nVerts = 4;
        face[3].vert = new VertexID[face[3].nVerts];
        face[3].vert[0].vertIndex = 8;
        face[3].vert[1].vertIndex = 15;
        face[3].vert[2].vertIndex = 12;
        face[3].vert[3].vertIndex = 11;

        // left bar
        // front
        face[4].nVerts = 4;
        face[4].vert = new VertexID[face[4].nVerts];
        face[4].vert[0].vertIndex = 0;
        face[4].vert[1].vertIndex = 1;
        face[4].vert[2].vertIndex = 2;
        face[4].vert[3].vertIndex = 3;
        // back
        face[5].nVerts = 4;
        face[5].vert = new VertexID[face[5].nVerts];
        face[5].vert[0].vertIndex = 11;
        face[5].vert[1].vertIndex = 10;
        face[5].vert[2].vertIndex = 9;
        face[5].vert[3].vertIndex = 8;
        // left
        face[6].nVerts = 4;
        face[6].vert = new VertexID[face[6].nVerts];
        face[6].vert[0].vertIndex = 9;
        face[6].vert[1].vertIndex = 1;
        face[6].vert[2].vertIndex = 0;
        face[6].vert[3].vertIndex = 8;
        // right
        face[7].nVerts = 4;
        face[7].vert = new VertexID[face[7].nVerts];
        face[7].vert[0].vertIndex = 3;
        face[7].vert[1].vertIndex = 2;
        face[7].vert[2].vertIndex = 10;
        face[7].vert[3].vertIndex = 11;
        // top
        face[12].nVerts = 4;
        face[12].vert = new VertexID[face[12].nVerts];
        face[12].vert[0].vertIndex = 1;
        face[12].vert[1].vertIndex = 9;
        face[12].vert[2].vertIndex = 10;
        face[12].vert[3].vertIndex = 2;        

        // right bar
        // front
        face[8].nVerts = 4;
        face[8].vert = new VertexID[face[8].nVerts];
        face[8].vert[0].vertIndex = 4;
        face[8].vert[1].vertIndex = 5;
        face[8].vert[2].vertIndex = 6;
        face[8].vert[3].vertIndex = 7;
        // back
        face[9].nVerts = 4;
        face[9].vert = new VertexID[face[9].nVerts];
        face[9].vert[0].vertIndex = 14;
        face[9].vert[1].vertIndex = 13;
        face[9].vert[2].vertIndex = 12;
        face[9].vert[3].vertIndex = 15;
        // left
        face[10].nVerts = 4;
        face[10].vert = new VertexID[face[10].nVerts];
        face[10].vert[0].vertIndex = 5;
        face[10].vert[1].vertIndex = 4;
        face[10].vert[2].vertIndex = 12;
        face[10].vert[3].vertIndex = 13;
        // right
        face[11].nVerts = 4;
        face[11].vert = new VertexID[face[11].nVerts];
        face[11].vert[0].vertIndex = 7;
        face[11].vert[1].vertIndex = 6;
        face[11].vert[2].vertIndex = 14;
        face[11].vert[3].vertIndex = 15;
        // top
        face[13].nVerts = 4;
        face[13].vert = new VertexID[face[13].nVerts];
        face[13].vert[0].vertIndex = 5;
        face[13].vert[1].vertIndex = 13;
        face[13].vert[2].vertIndex = 14;
        face[13].vert[3].vertIndex = 6; 


        CalculateFacesNorm();
    }
};

class Nut: public Mesh {
public:
    float size, height, radius;
    int nSegments;
    Nut(float size, float height, float radius, int nSegments): Mesh(), 
        size(size), height(height), radius(radius), nSegments(nSegments) {}
    ~Nut() {}

    virtual void init() {
        // init vertices
        numVerts = 16 + 4*nSegments;
        pt = new Point3f[numVerts];

        // left rectangular
        pt[numVerts - 1].set(-size/2, -size/2, height);
        pt[numVerts - 5].set(-radius, -size/2, height);
        pt[numVerts - 8].set(-radius, size/2, height);
        pt[numVerts - 4].set(-size/2, size/2, height);

        pt[numVerts - 9].set(-size/2, -size/2, 0);
        pt[numVerts - 13].set(-radius, -size/2, 0);
        pt[numVerts - 16].set(-radius, size/2, 0);
        pt[numVerts - 12].set(-size/2, size/2, 0);

        // right rectangular
        pt[numVerts - 6].set(radius, -size/2, height);
        pt[numVerts - 2].set(size/2, -size/2, height);
        pt[numVerts - 3].set(size/2, size/2, height);
        pt[numVerts - 7].set(radius, size/2, height);

        pt[numVerts - 14].set(radius, -size/2, 0);
        pt[numVerts - 10].set(size/2, -size/2, 0);
        pt[numVerts - 11].set(size/2, size/2, 0);
        pt[numVerts - 15].set(radius, size/2, 0);

        // circular vertices
        float angle = 0;
        float deltaAngle = 2*PI / nSegments;
        for(int i=0; i<nSegments; i++) {
            // front
            pt[i].set(radius*cos(angle), radius*sin(angle), height);
            // back
            pt[i + nSegments].set(radius*cos(angle), radius*sin(angle), 0);

            if(i < nSegments/2) {
                // front
                pt[i + 2*nSegments].set(radius*cos(angle), size/2, height);
                // back
                pt[i + 2*nSegments + nSegments/2].set(radius*cos(angle), size/2, 0);
            } else {
                // front
                pt[i + 3*nSegments - nSegments/2].set(radius*cos(angle), -size/2, height);
                // back
                pt[i + 3*nSegments].set(radius*cos(angle), -size/2, 0);             
            }
            angle += deltaAngle;
        }


        // init faces
        numFaces = 8 + 3*nSegments;
        face = new Face[numFaces];

        // left rectangular
        // front
        face[numFaces - 1].nVerts = 4;
        face[numFaces - 1].vert = new VertexID[face[numFaces - 1].nVerts];
        face[numFaces - 1].vert[0].vertIndex = numVerts - 1;
        face[numFaces - 1].vert[1].vertIndex = numVerts - 4;
        face[numFaces - 1].vert[2].vertIndex = numVerts - 8;
        face[numFaces - 1].vert[3].vertIndex = numVerts - 5;
        // back
        face[numFaces - 2].nVerts = 4;
        face[numFaces - 2].vert = new VertexID[face[numFaces - 2].nVerts];
        face[numFaces - 2].vert[0].vertIndex = numVerts - 9;
        face[numFaces - 2].vert[1].vertIndex = numVerts - 13;
        face[numFaces - 2].vert[2].vertIndex = numVerts - 16;
        face[numFaces - 2].vert[3].vertIndex = numVerts - 12;
        // left
        face[numFaces - 3].nVerts = 4;
        face[numFaces - 3].vert = new VertexID[face[numFaces - 3].nVerts];
        face[numFaces - 3].vert[0].vertIndex = numVerts - 4;
        face[numFaces - 3].vert[1].vertIndex = numVerts - 1;
        face[numFaces - 3].vert[2].vertIndex = numVerts - 9;
        face[numFaces - 3].vert[3].vertIndex = numVerts - 12;

        // right rectangular
        // front
        face[numFaces - 4].nVerts = 4;
        face[numFaces - 4].vert = new VertexID[face[numFaces - 4].nVerts];
        face[numFaces - 4].vert[0].vertIndex = numVerts - 2;
        face[numFaces - 4].vert[1].vertIndex = numVerts - 6;
        face[numFaces - 4].vert[2].vertIndex = numVerts - 7;
        face[numFaces - 4].vert[3].vertIndex = numVerts - 3;
        // back
        face[numFaces - 5].nVerts = 4;
        face[numFaces - 5].vert = new VertexID[face[numFaces - 5].nVerts];
        face[numFaces - 5].vert[0].vertIndex = numVerts - 14;
        face[numFaces - 5].vert[1].vertIndex = numVerts - 10;
        face[numFaces - 5].vert[2].vertIndex = numVerts - 11;
        face[numFaces - 5].vert[3].vertIndex = numVerts - 15;
        // right
        face[numFaces - 6].nVerts = 4;
        face[numFaces - 6].vert = new VertexID[face[numFaces - 6].nVerts];
        face[numFaces - 6].vert[0].vertIndex = numVerts - 2;
        face[numFaces - 6].vert[1].vertIndex = numVerts - 3;
        face[numFaces - 6].vert[2].vertIndex = numVerts - 11;
        face[numFaces - 6].vert[3].vertIndex = numVerts - 10;


        // top
        face[numFaces - 7].nVerts = 4;
        face[numFaces - 7].vert = new VertexID[face[numFaces - 7].nVerts];
        face[numFaces - 7].vert[0].vertIndex = numVerts - 3;
        face[numFaces - 7].vert[1].vertIndex = numVerts - 4;
        face[numFaces - 7].vert[2].vertIndex = numVerts - 12;
        face[numFaces - 7].vert[3].vertIndex = numVerts - 11;
        // bottom
        face[numFaces - 8].nVerts = 4;
        face[numFaces - 8].vert = new VertexID[face[numFaces - 8].nVerts];
        face[numFaces - 8].vert[0].vertIndex = numVerts - 1;
        face[numFaces - 8].vert[1].vertIndex = numVerts - 2;
        face[numFaces - 8].vert[2].vertIndex = numVerts - 10;
        face[numFaces - 8].vert[3].vertIndex = numVerts - 9;

        
        // circular faces
        for(int i=0; i<nSegments; i++) {
            face[i].nVerts = 4;
            face[i].vert = new VertexID[face[i].nVerts];
            face[i].vert[0].vertIndex = i % nSegments;
            face[i].vert[1].vertIndex = (i+1) % nSegments;
            face[i].vert[2].vertIndex = (i+1) % nSegments + nSegments;
            face[i].vert[3].vertIndex = i % nSegments + nSegments;

            if(i < nSegments/2) {
                bool isLast = (i == nSegments/2 - 1);
                // front
                face[i + nSegments].nVerts = 4;
                face[i + nSegments].vert = new VertexID[face[i + nSegments].nVerts];
                face[i + nSegments].vert[0].vertIndex = i % nSegments;
                face[i + nSegments].vert[1].vertIndex = (i+1) % nSegments;
                face[i + nSegments].vert[2].vertIndex = isLast ? numVerts - 8 : (i+1) % nSegments + 2*nSegments;
                face[i + nSegments].vert[3].vertIndex = i % nSegments + 2*nSegments;
                // back
                face[i + nSegments + nSegments/2].nVerts = 4;
                face[i + nSegments + nSegments/2].vert = new VertexID[face[i + nSegments + nSegments/2].nVerts];
                face[i + nSegments + nSegments/2].vert[0].vertIndex = (i+1) % nSegments + nSegments;
                face[i + nSegments + nSegments/2].vert[1].vertIndex = i % nSegments + nSegments;
                face[i + nSegments + nSegments/2].vert[2].vertIndex = i % nSegments + 2*nSegments + nSegments/2;
                face[i + nSegments + nSegments/2].vert[3].vertIndex = isLast ? numVerts - 16 : (i+1) % nSegments + 2*nSegments + nSegments/2;
            } else {
                bool isLast = (i == nSegments - 1);
                // front
                face[i + 2*nSegments - nSegments/2].nVerts = 4;
                face[i + 2*nSegments - nSegments/2].vert = new VertexID[face[i + 2*nSegments - nSegments/2].nVerts];
                face[i + 2*nSegments - nSegments/2].vert[0].vertIndex = i % nSegments;
                face[i + 2*nSegments - nSegments/2].vert[1].vertIndex = (i+1) % nSegments;
                face[i + 2*nSegments - nSegments/2].vert[2].vertIndex = isLast ? numVerts - 6 : (i+1) % nSegments + 3*nSegments - nSegments/2;
                face[i + 2*nSegments - nSegments/2].vert[3].vertIndex = i % nSegments + 3*nSegments - nSegments/2;
                // back
                face[i + 2*nSegments].nVerts = 4;
                face[i + 2*nSegments].vert = new VertexID[face[i + 2*nSegments].nVerts];
                face[i + 2*nSegments].vert[0].vertIndex = (i+1) % nSegments + nSegments;
                face[i + 2*nSegments].vert[1].vertIndex = i % nSegments + nSegments;
                face[i + 2*nSegments].vert[2].vertIndex = i % nSegments + 3*nSegments;
                face[i + 2*nSegments].vert[3].vertIndex = isLast ? numVerts - 14 : (i+1) % nSegments + 3*nSegments;              
            }
        }

        CalculateFacesNorm();
    }
};

class LinkBar: public Mesh {
public:
    float length, width, height, slotLength, slotShift, slotWidth;
    int nSegments;
    LinkBar(float length, float width, float height, float slotLength, float slotShift, float slotWidth, int nSegments):
        Mesh(), length(length), width(width), height(height),
        slotLength(slotLength), slotShift(slotShift), slotWidth(slotWidth), nSegments(nSegments) {}
    ~LinkBar() {}

    virtual void init() {
        float radius = width/2;
        // init vertices
        numVerts = 28 + 4*nSegments;
        pt = new Point3f[numVerts];
        
        // front
        pt[numVerts - 1].set(-length/2 + radius, 0, height);
        pt[numVerts - 2].set(-length/2 + radius, -width/2, height);
        pt[numVerts - 3].set(-length/2 + radius, width/2, height);
        pt[numVerts - 4].set(-slotLength/2 + slotShift, width/2, height);
        pt[numVerts - 5].set(slotLength/2 + slotShift, width/2, height);
        pt[numVerts - 6].set(length/2 - radius, width/2, height);
        pt[numVerts - 7].set(length/2 - radius, -width/2, height);
        pt[numVerts - 8].set(slotLength/2 + slotShift, -width/2, height);
        pt[numVerts - 9].set(-slotLength/2 + slotShift, -width/2, height);
        pt[numVerts - 10].set(-slotLength/2 + slotShift, slotWidth/2, height);
        pt[numVerts - 11].set(-slotLength/2 + slotShift, -slotWidth/2, height);
        pt[numVerts - 12].set(slotLength/2 + slotShift, slotWidth/2, height);
        pt[numVerts - 13].set(slotLength/2 + slotShift, -slotWidth/2, height);
        pt[numVerts - 14].set(length/2 - radius, 0, height);

        // back
        pt[numVerts - 15].set(-length/2 + radius, 0, 0);
        pt[numVerts - 16].set(-length/2 + radius, -width/2, 0);
        pt[numVerts - 17].set(-length/2 + radius, width/2, 0);
        pt[numVerts - 18].set(-slotLength/2 + slotShift, width/2, 0);
        pt[numVerts - 19].set(slotLength/2 + slotShift, width/2, 0);
        pt[numVerts - 20].set(length/2 - radius, width/2, 0);
        pt[numVerts - 21].set(length/2 - radius, -width/2, 0);
        pt[numVerts - 22].set(slotLength/2 + slotShift, -width/2, 0);
        pt[numVerts - 23].set(-slotLength/2 + slotShift, -width/2, 0);
        pt[numVerts - 24].set(-slotLength/2 + slotShift, slotWidth/2, 0);
        pt[numVerts - 25].set(-slotLength/2 + slotShift, -slotWidth/2, 0);
        pt[numVerts - 26].set(slotLength/2 + slotShift, slotWidth/2, 0);
        pt[numVerts - 27].set(slotLength/2 + slotShift, -slotWidth/2, 0);
        pt[numVerts - 28].set(length/2 - radius, 0, 0);

    
        float angle = 0;
        float deltaAngle = PI / nSegments;
        for(int i=0; i<nSegments; i++) {
            // left circular
            // front
            pt[i].set(-length/2 + radius - radius*sin(angle), radius*cos(angle), height);
            // back
            pt[i+nSegments].set(-length/2 + radius - radius*sin(angle), radius*cos(angle), 0);

            // right circular
            // front
            pt[i+2*nSegments].set(length/2 - radius + radius*sin(angle), radius*cos(angle), height);
            // back
            pt[i+3*nSegments].set(length/2 - radius + radius*sin(angle), radius*cos(angle), 0);

            angle += deltaAngle;
        }


        // init faces
        numFaces = 18 + 6*nSegments;
        face = new Face[numFaces];

        // left rectangular
        // front
        face[numFaces-1].nVerts = 4;
        face[numFaces-1].vert = new VertexID[face[numFaces-1].nVerts];
        face[numFaces-1].vert[0].vertIndex = numVerts-2;
        face[numFaces-1].vert[1].vertIndex = numVerts-3;
        face[numFaces-1].vert[2].vertIndex = numVerts-4;
        face[numFaces-1].vert[3].vertIndex = numVerts-9;
        // back
        face[numFaces-2].nVerts = 4;
        face[numFaces-2].vert = new VertexID[face[numFaces-2].nVerts];
        face[numFaces-2].vert[0].vertIndex = numVerts-18;
        face[numFaces-2].vert[1].vertIndex = numVerts-17;
        face[numFaces-2].vert[2].vertIndex = numVerts-16;
        face[numFaces-2].vert[3].vertIndex = numVerts-23;
        // top
        face[numFaces-3].nVerts = 4;
        face[numFaces-3].vert = new VertexID[face[numFaces-3].nVerts];
        face[numFaces-3].vert[0].vertIndex = numVerts-3;
        face[numFaces-3].vert[1].vertIndex = numVerts-17;
        face[numFaces-3].vert[2].vertIndex = numVerts-18;
        face[numFaces-3].vert[3].vertIndex = numVerts-4;
        // bottom
        face[numFaces-4].nVerts = 4;
        face[numFaces-4].vert = new VertexID[face[numFaces-4].nVerts];
        face[numFaces-4].vert[0].vertIndex = numVerts-2;
        face[numFaces-4].vert[1].vertIndex = numVerts-9;
        face[numFaces-4].vert[2].vertIndex = numVerts-23;
        face[numFaces-4].vert[3].vertIndex = numVerts-16;
        // right
        face[numFaces-17].nVerts = 4;
        face[numFaces-17].vert = new VertexID[face[numFaces-17].nVerts];
        face[numFaces-17].vert[0].vertIndex = numVerts-9;
        face[numFaces-17].vert[1].vertIndex = numVerts-4;
        face[numFaces-17].vert[2].vertIndex = numVerts-18;
        face[numFaces-17].vert[3].vertIndex = numVerts-23;


        // top rectangular
        // front
        face[numFaces-5].nVerts = 4;
        face[numFaces-5].vert = new VertexID[face[numFaces-5].nVerts];
        face[numFaces-5].vert[0].vertIndex = numVerts-10;
        face[numFaces-5].vert[1].vertIndex = numVerts-4;
        face[numFaces-5].vert[2].vertIndex = numVerts-5;
        face[numFaces-5].vert[3].vertIndex = numVerts-12;
        // back
        face[numFaces-6].nVerts = 4;
        face[numFaces-6].vert = new VertexID[face[numFaces-6].nVerts];
        face[numFaces-6].vert[0].vertIndex = numVerts-18;
        face[numFaces-6].vert[1].vertIndex = numVerts-24;
        face[numFaces-6].vert[2].vertIndex = numVerts-26;
        face[numFaces-6].vert[3].vertIndex = numVerts-19;
        // top
        face[numFaces-7].nVerts = 4;
        face[numFaces-7].vert = new VertexID[face[numFaces-7].nVerts];
        face[numFaces-7].vert[0].vertIndex = numVerts-4;
        face[numFaces-7].vert[1].vertIndex = numVerts-18;
        face[numFaces-7].vert[2].vertIndex = numVerts-19;
        face[numFaces-7].vert[3].vertIndex = numVerts-5;
        // bottom
        face[numFaces-8].nVerts = 4;
        face[numFaces-8].vert = new VertexID[face[numFaces-8].nVerts];
        face[numFaces-8].vert[0].vertIndex = numVerts-10;
        face[numFaces-8].vert[1].vertIndex = numVerts-12;
        face[numFaces-8].vert[2].vertIndex = numVerts-26;
        face[numFaces-8].vert[3].vertIndex = numVerts-24;


        // bottom rectangular
        // front
        face[numFaces-9].nVerts = 4;
        face[numFaces-9].vert = new VertexID[face[numFaces-9].nVerts];
        face[numFaces-9].vert[0].vertIndex = numVerts-9;
        face[numFaces-9].vert[1].vertIndex = numVerts-11;
        face[numFaces-9].vert[2].vertIndex = numVerts-13;
        face[numFaces-9].vert[3].vertIndex = numVerts-8;
        // back
        face[numFaces-10].nVerts = 4;
        face[numFaces-10].vert = new VertexID[face[numFaces-10].nVerts];
        face[numFaces-10].vert[0].vertIndex = numVerts-27;
        face[numFaces-10].vert[1].vertIndex = numVerts-25;
        face[numFaces-10].vert[2].vertIndex = numVerts-23;
        face[numFaces-10].vert[3].vertIndex = numVerts-22;
        // top
        face[numFaces-11].nVerts = 4;
        face[numFaces-11].vert = new VertexID[face[numFaces-11].nVerts];
        face[numFaces-11].vert[0].vertIndex = numVerts-11;
        face[numFaces-11].vert[1].vertIndex = numVerts-25;
        face[numFaces-11].vert[2].vertIndex = numVerts-27;
        face[numFaces-11].vert[3].vertIndex = numVerts-13;
        // bottom
        face[numFaces-12].nVerts = 4;
        face[numFaces-12].vert = new VertexID[face[numFaces-12].nVerts];
        face[numFaces-12].vert[0].vertIndex = numVerts-9;
        face[numFaces-12].vert[1].vertIndex = numVerts-8;
        face[numFaces-12].vert[2].vertIndex = numVerts-22;
        face[numFaces-12].vert[3].vertIndex = numVerts-23;


        // right rectangular
        // front
        face[numFaces-13].nVerts = 4;
        face[numFaces-13].vert = new VertexID[face[numFaces-13].nVerts];
        face[numFaces-13].vert[0].vertIndex = numVerts-5;
        face[numFaces-13].vert[1].vertIndex = numVerts-6;
        face[numFaces-13].vert[2].vertIndex = numVerts-7;
        face[numFaces-13].vert[3].vertIndex = numVerts-8;
        // back
        face[numFaces-14].nVerts = 4;
        face[numFaces-14].vert = new VertexID[face[numFaces-14].nVerts];
        face[numFaces-14].vert[0].vertIndex = numVerts-19;
        face[numFaces-14].vert[1].vertIndex = numVerts-22;
        face[numFaces-14].vert[2].vertIndex = numVerts-21;
        face[numFaces-14].vert[3].vertIndex = numVerts-20;
        // top
        face[numFaces-15].nVerts = 4;
        face[numFaces-15].vert = new VertexID[face[numFaces-15].nVerts];
        face[numFaces-15].vert[0].vertIndex = numVerts-6;
        face[numFaces-15].vert[1].vertIndex = numVerts-5;
        face[numFaces-15].vert[2].vertIndex = numVerts-19;
        face[numFaces-15].vert[3].vertIndex = numVerts-20;
        // bottom
        face[numFaces-16].nVerts = 4;
        face[numFaces-16].vert = new VertexID[face[numFaces-16].nVerts];
        face[numFaces-16].vert[0].vertIndex = numVerts-8;
        face[numFaces-16].vert[1].vertIndex = numVerts-7;
        face[numFaces-16].vert[2].vertIndex = numVerts-21;
        face[numFaces-16].vert[3].vertIndex = numVerts-22;
        // left
        face[numFaces-18].nVerts = 4;
        face[numFaces-18].vert = new VertexID[face[numFaces-18].nVerts];
        face[numFaces-18].vert[0].vertIndex = numVerts-5;
        face[numFaces-18].vert[1].vertIndex = numVerts-8;
        face[numFaces-18].vert[2].vertIndex = numVerts-22;
        face[numFaces-18].vert[3].vertIndex = numVerts-19;


        for(int i=0; i<nSegments; i++) {
            // left circular
            // front
            face[i].nVerts = 3;
            face[i].vert = new VertexID[face[i].nVerts];
            face[i].vert[0].vertIndex = numVerts-1;
            face[i].vert[1].vertIndex = (i == nSegments-1) ?  numVerts-2 : i+1;
            face[i].vert[2].vertIndex = i;
            // back
            face[i+nSegments].nVerts = 3;
            face[i+nSegments].vert = new VertexID[face[i+nSegments].nVerts];
            face[i+nSegments].vert[0].vertIndex = numVerts-15;
            face[i+nSegments].vert[1].vertIndex = i + nSegments;
            face[i+nSegments].vert[2].vertIndex = (i == nSegments-1) ?  numVerts-16 : i+1 + nSegments;
            // connect faces
            face[i+2*nSegments].nVerts = 4;
            face[i+2*nSegments].vert = new VertexID[face[i+2*nSegments].nVerts];
            face[i+2*nSegments].vert[0].vertIndex = i;
            face[i+2*nSegments].vert[1].vertIndex = (i == nSegments-1) ? numVerts-2 : i+1;
            face[i+2*nSegments].vert[2].vertIndex = (i == nSegments-1) ? numVerts-16 : i+1 + nSegments;
            face[i+2*nSegments].vert[3].vertIndex = i + nSegments;     

            // right circular
            // front
            face[i+3*nSegments].nVerts = 3;
            face[i+3*nSegments].vert = new VertexID[face[i+3*nSegments].nVerts];
            face[i+3*nSegments].vert[0].vertIndex = numVerts-14;
            face[i+3*nSegments].vert[1].vertIndex = i + 2*nSegments;
            face[i+3*nSegments].vert[2].vertIndex = (i == nSegments-1) ?  numVerts-7 : i+1 + 2*nSegments;
            // back
            face[i+4*nSegments].nVerts = 3;
            face[i+4*nSegments].vert = new VertexID[face[i+4*nSegments].nVerts];
            face[i+4*nSegments].vert[0].vertIndex = numVerts-28;
            face[i+4*nSegments].vert[1].vertIndex = (i == nSegments-1) ?  numVerts-21 : i+1 + 3*nSegments;
            face[i+4*nSegments].vert[2].vertIndex = i + 3*nSegments;
            // connect faces
            face[i+5*nSegments].nVerts = 4;
            face[i+5*nSegments].vert = new VertexID[face[i+5*nSegments].nVerts];
            face[i+5*nSegments].vert[0].vertIndex = (i == nSegments-1) ? numVerts-7 : i+1 + 2*nSegments;
            face[i+5*nSegments].vert[1].vertIndex = i + 2*nSegments;
            face[i+5*nSegments].vert[2].vertIndex = i + 3*nSegments; 
            face[i+5*nSegments].vert[3].vertIndex = (i == nSegments-1) ? numVerts-21 : i+1 + 3*nSegments;
        }


        CalculateFacesNorm();
    }
};

class FirstTile: public Tile {
public:
    float size;
    int nSegments;
    FirstTile(float size, int nSegments): Tile(), size(size), nSegments(nSegments) {}
    ~FirstTile() {}

    virtual void init() {
        float radius = 0.6*size;
        float innerradius = 0.7*radius;
        // init vertices
        numVerts = 4 + 8*(nSegments+1);
        pt = new Point3f[numVerts];

        // square vertices
        pt[numVerts - 1].set(-size/2, size/2, 0);
        pt[numVerts - 2].set(size/2, size/2, 0);
        pt[numVerts - 3].set(size/2, -size/2, 0);
        pt[numVerts - 4].set(-size/2, -size/2, 0);

        // circular vertices
        float angle = 0;
        float deltaAngle = PI / (2*nSegments);
        for(int i=0; i<nSegments+1; i++) {
            //front
            // bottom left vertices
            // inner circle
            pt[i].set(-size/2 + innerradius*cos(angle), -size/2 + innerradius*sin(angle), 0.001);
            //outer circle
            pt[i + nSegments+1].set(-size/2 + radius*cos(angle), -size/2 + radius*sin(angle), 0.001);

            // top right vertices
            // inner circle
            pt[i + 2*(nSegments+1)].set(size/2 - innerradius*cos(angle), size/2 - innerradius*sin(angle), 0.001);
            //outer circle
            pt[i + 3*(nSegments+1)].set(size/2 - radius*cos(angle), size/2 - radius*sin(angle), 0.001);

            // back
            // bottom left vertices
            // inner circle
            pt[i + 4*(nSegments+1)].set(-size/2 + innerradius*cos(angle), -size/2 + innerradius*sin(angle), -0.001);
            //outer circle
            pt[i + 5*(nSegments+1)].set(-size/2 + radius*cos(angle), -size/2 + radius*sin(angle), -0.001);

            // top right vertices
            // inner circle
            pt[i + 6*(nSegments+1)].set(size/2 - innerradius*cos(angle), size/2 - innerradius*sin(angle), -0.001);
            //outer circle
            pt[i + 7*(nSegments+1)].set(size/2 - radius*cos(angle), size/2 - radius*sin(angle), -0.001);

            angle += deltaAngle;
        }


        // init faces
        numFaces = 1 + 4*nSegments;
        face = new Face[numFaces];

        // square face
        face[numFaces - 1].nVerts = 4;
        face[numFaces - 1].vert = new VertexID[face[numFaces - 1].nVerts];
        face[numFaces - 1].vert[0].vertIndex = numVerts - 1;
        face[numFaces - 1].vert[1].vertIndex = numVerts - 2;
        face[numFaces - 1].vert[2].vertIndex = numVerts - 3;
        face[numFaces - 1].vert[3].vertIndex = numVerts - 4;

        // circular faces
        for(int i=0; i<nSegments; i++) {
            // front
            // bottom left
            face[i].nVerts = 4;
            face[i].vert = new VertexID[face[i].nVerts];
            face[i].vert[0].vertIndex = i;
            face[i].vert[1].vertIndex = i+1;
            face[i].vert[2].vertIndex = i+1 + nSegments+1;
            face[i].vert[3].vertIndex = i + nSegments+1;

            // top right
            face[i + nSegments].nVerts = 4;
            face[i + nSegments].vert = new VertexID[face[i + nSegments].nVerts];
            face[i + nSegments].vert[0].vertIndex = i + 3*(nSegments+1);
            face[i + nSegments].vert[1].vertIndex = i + 2*(nSegments+1);
            face[i + nSegments].vert[2].vertIndex = i+1 + 2*(nSegments+1);
            face[i + nSegments].vert[3].vertIndex = i+1 + 3*(nSegments+1);

            // back
            // bottom left
            face[i + 2*nSegments].nVerts = 4;
            face[i + 2*nSegments].vert = new VertexID[face[i + 2*nSegments].nVerts];
            face[i + 2*nSegments].vert[0].vertIndex = i + 4*(nSegments+1);
            face[i + 2*nSegments].vert[1].vertIndex = i + 5*(nSegments+1);
            face[i + 2*nSegments].vert[2].vertIndex = i+1 + 5*(nSegments+1);
            face[i + 2*nSegments].vert[3].vertIndex = i+1 + 4*(nSegments+1);

            // top right
            face[i + 3*nSegments].nVerts = 4;
            face[i + 3*nSegments].vert = new VertexID[face[i + 3*nSegments].nVerts];
            face[i + 3*nSegments].vert[0].vertIndex = i + 6*(nSegments+1);
            face[i + 3*nSegments].vert[1].vertIndex = i + 7*(nSegments+1);
            face[i + 3*nSegments].vert[2].vertIndex = i+1 + 7*(nSegments+1);
            face[i + 3*nSegments].vert[3].vertIndex = i+1 + 6*(nSegments+1);
        }
    }

};

class SecondTile: public Tile {
public:
    float size;
    SecondTile(float size): Tile(), size(size) {}
    ~SecondTile() {}

    virtual void init() {
        float patternSize = size/5;
        // init vertices
        numVerts = 16;
        pt = new Point3f[numVerts];
        pt[0].set(-size/2, size/2, 0);
        pt[1].set(size/2, size/2, 0);
        pt[2].set(-size/2, -size/2, 0);
        pt[3].set(size/2, -size/2, 0);
        pt[4].set(-patternSize/2, size/2, 0);
        pt[5].set(patternSize/2, size/2, 0);
        pt[6].set(size/2, patternSize/2, 0);
        pt[7].set(size/2, -patternSize/2, 0);
        pt[8].set(patternSize/2, -size/2, 0);
        pt[9].set(-patternSize/2, -size/2, 0);
        pt[10].set(-size/2, -patternSize/2, 0);
        pt[11].set(-size/2, patternSize/2, 0);
        pt[12].set(-patternSize/2, patternSize/2, 0);
        pt[13].set(patternSize/2, patternSize/2, 0);
        pt[14].set(patternSize/2, -patternSize/2, 0);
        pt[15].set(-patternSize/2, -patternSize/2, 0);

        // init faces
        numFaces = 6;
        face = new Face[numFaces];

        // top left
        face[0].nVerts = 4;
        face[0].vert = new VertexID[face[0].nVerts];
        face[0].vert[0].vertIndex = 0;
        face[0].vert[1].vertIndex = 4;
        face[0].vert[2].vertIndex = 12;
        face[0].vert[3].vertIndex = 11;

        // top right
        face[1].nVerts = 4;
        face[1].vert = new VertexID[face[1].nVerts];
        face[1].vert[0].vertIndex = 5;
        face[1].vert[1].vertIndex = 1;
        face[1].vert[2].vertIndex = 6;
        face[1].vert[3].vertIndex = 13;

        // bottom left
        face[2].nVerts = 4;
        face[2].vert = new VertexID[face[2].nVerts];
        face[2].vert[0].vertIndex = 2;
        face[2].vert[1].vertIndex = 10;
        face[2].vert[2].vertIndex = 15;
        face[2].vert[3].vertIndex = 9;

        // bottom right
        face[3].nVerts = 4;
        face[3].vert = new VertexID[face[3].nVerts];
        face[3].vert[0].vertIndex = 8;
        face[3].vert[1].vertIndex = 14;
        face[3].vert[2].vertIndex = 7;
        face[3].vert[3].vertIndex = 3;

        // vertical
        face[4].nVerts = 4;
        face[4].vert = new VertexID[face[4].nVerts];
        face[4].vert[0].vertIndex = 4;
        face[4].vert[1].vertIndex = 5;
        face[4].vert[2].vertIndex = 8;
        face[4].vert[3].vertIndex = 9;

        // horizontal
        face[5].nVerts = 4;
        face[5].vert = new VertexID[face[5].nVerts];
        face[5].vert[0].vertIndex = 11;
        face[5].vert[1].vertIndex = 6;
        face[5].vert[2].vertIndex = 7;
        face[5].vert[3].vertIndex = 10;
    }

    void DrawColor() {
        glDisable(GL_LIGHTING);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        for (int f = 0; f < numFaces; f++)
        {
            if (f == numFaces-1 || f == numFaces-2) {
                glColor3f(secondColor.red, secondColor.green, secondColor.blue);
            } else {
                glColor3f(firstColor.red, firstColor.green, firstColor.blue);
            }
            
            glBegin(GL_POLYGON);
            for (int v = 0; v < face[f].nVerts; v++)
            {
                int	iv = face[f].vert[v].vertIndex;
                // int	ic = face[f].vert[v].colorIndex; 
                glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
            }
            glEnd();
        }

        glEnable(GL_LIGHTING);
    }
};

class ThirdTile: public Tile {
public:
    float size;
    int nSegments;
    ThirdTile(float size, int nSegments): Tile(), size(size), nSegments(nSegments) {}
    ~ThirdTile() {}

    virtual void init() {
        float radius = 0.8*size/2;
        float patternSize = size/5;
        float innerradius = radius - patternSize;
        // init vertices
        numVerts = 4 + 2*(8 + 2*nSegments + 8);
        pt = new Point3f[numVerts];

        // square
        pt[numVerts-1].set(-size/2, size/2, 0);
        pt[numVerts-2].set(size/2, size/2, 0);
        pt[numVerts-3].set(size/2, -size/2, 0);
        pt[numVerts-4].set(-size/2, -size/2, 0);



        // vertical
        // front
        // top
        pt[numVerts - 5].set(-patternSize/2, size/2, 0.001);
        pt[numVerts - 6].set(patternSize/2, size/2, 0.001);
        pt[numVerts - 21].set(-patternSize/2, innerradius+0.001, 0.001);
        pt[numVerts - 22].set(patternSize/2, innerradius+0.001, 0.001);
        // bottom
        pt[numVerts - 7].set(patternSize/2, -size/2, 0.001);
        pt[numVerts - 8].set(-patternSize/2, -size/2, 0.001);
        pt[numVerts - 23].set(-patternSize/2, -(innerradius+0.001), 0.001);
        pt[numVerts - 24].set(patternSize/2, -(innerradius+0.001), 0.001);        

        // back
        // top
        pt[numVerts - 13].set(-patternSize/2, size/2, -0.001);
        pt[numVerts - 14].set(patternSize/2, size/2, -0.001);
        pt[numVerts - 25].set(-patternSize/2, innerradius+0.001, -0.001);
        pt[numVerts - 26].set(patternSize/2, innerradius+0.001, -0.001);
        // bottom
        pt[numVerts - 15].set(patternSize/2, -size/2, -0.001);
        pt[numVerts - 16].set(-patternSize/2, -size/2, -0.001);
        pt[numVerts - 27].set(-patternSize/2, -(innerradius+0.001), -0.001);
        pt[numVerts - 28].set(patternSize/2, -(innerradius+0.001), -0.001);         




        // horizontal
        // front
        // left
        pt[numVerts - 9].set(-size/2, patternSize/2, 0.001);
        pt[numVerts - 10].set(-size/2, -patternSize/2, 0.001);
        pt[numVerts - 29].set(-(innerradius+0.001), patternSize/2, 0.001);
        pt[numVerts - 30].set(-(innerradius+0.001), -patternSize/2, 0.001);
        // right
        pt[numVerts - 11].set(size/2, patternSize/2, 0.001);
        pt[numVerts - 12].set(size/2, -patternSize/2, 0.001);
        pt[numVerts - 31].set((innerradius+0.001), patternSize/2, 0.001);
        pt[numVerts - 32].set((innerradius+0.001), -patternSize/2, 0.001);

        // back
        // left
        pt[numVerts - 17].set(-size/2, patternSize/2, -0.001);
        pt[numVerts - 18].set(-size/2, -patternSize/2, -0.001);
        pt[numVerts - 33].set(-(innerradius+0.001), patternSize/2, -0.001);
        pt[numVerts - 34].set(-(innerradius+0.001), -patternSize/2, -0.001);
        // right
        pt[numVerts - 19].set(size/2, patternSize/2, -0.001);
        pt[numVerts - 20].set(size/2, -patternSize/2, -0.001);
        pt[numVerts - 35].set((innerradius+0.001), patternSize/2, -0.001);
        pt[numVerts - 36].set((innerradius+0.001), -patternSize/2, -0.001);




        // circular faces
        float angle = 0;
        float deltaAngle = 2*PI / nSegments;
        for(int i=0; i<nSegments; i++) {
            // outer
            // front
            pt[i].set(radius*cos(angle), radius*sin(angle), 0.001);
            // back
            pt[i + nSegments].set(radius*cos(angle), radius*sin(angle), -0.001);
            
            // inner
            // front
            pt[i + 2*nSegments].set(innerradius*cos(angle), innerradius*sin(angle), 0.001);
            // back
            pt[i + 3*nSegments].set(innerradius*cos(angle), innerradius*sin(angle), -0.001);
            angle += deltaAngle;
        }


        // init faces
        numFaces = 1 + 2*(4+nSegments);
        face = new Face[numFaces];

        // square
        face[numFaces - 1].nVerts = 4;
        face[numFaces - 1].vert = new VertexID[face[numFaces - 1].nVerts];
        face[numFaces - 1].vert[0].vertIndex = numVerts - 1;
        face[numFaces - 1].vert[1].vertIndex = numVerts - 2;
        face[numFaces - 1].vert[2].vertIndex = numVerts - 3;
        face[numFaces - 1].vert[3].vertIndex = numVerts - 4;

        // vertical
        // top
        // front
        face[numFaces - 2].nVerts = 4;
        face[numFaces - 2].vert = new VertexID[face[numFaces - 2].nVerts];
        face[numFaces - 2].vert[0].vertIndex = numVerts - 5;
        face[numFaces - 2].vert[1].vertIndex = numVerts - 6;
        face[numFaces - 2].vert[2].vertIndex = numVerts - 22;
        face[numFaces - 2].vert[3].vertIndex = numVerts - 21;
        // back
        face[numFaces - 3].nVerts = 4;
        face[numFaces - 3].vert = new VertexID[face[numFaces - 3].nVerts];
        face[numFaces - 3].vert[0].vertIndex = numVerts - 13;
        face[numFaces - 3].vert[1].vertIndex = numVerts - 25;
        face[numFaces - 3].vert[2].vertIndex = numVerts - 26;
        face[numFaces - 3].vert[3].vertIndex = numVerts - 14;
        
        // bottom
        // front
        face[numFaces - 4].nVerts = 4;
        face[numFaces - 4].vert = new VertexID[face[numFaces - 4].nVerts];
        face[numFaces - 4].vert[0].vertIndex = numVerts - 8;
        face[numFaces - 4].vert[1].vertIndex = numVerts - 23;
        face[numFaces - 4].vert[2].vertIndex = numVerts - 24;
        face[numFaces - 4].vert[3].vertIndex = numVerts - 7;
        // back
        face[numFaces - 5].nVerts = 4;
        face[numFaces - 5].vert = new VertexID[face[numFaces - 5].nVerts];
        face[numFaces - 5].vert[0].vertIndex = numVerts - 16;
        face[numFaces - 5].vert[1].vertIndex = numVerts - 15;
        face[numFaces - 5].vert[2].vertIndex = numVerts - 28;
        face[numFaces - 5].vert[3].vertIndex = numVerts - 27;



        // horizontal
        // left
        // front
        face[numFaces - 6].nVerts = 4;
        face[numFaces - 6].vert = new VertexID[face[numFaces - 6].nVerts];
        face[numFaces - 6].vert[0].vertIndex = numVerts - 9;
        face[numFaces - 6].vert[1].vertIndex = numVerts - 29;
        face[numFaces - 6].vert[2].vertIndex = numVerts - 30;
        face[numFaces - 6].vert[3].vertIndex = numVerts - 10;
        // back
        face[numFaces - 7].nVerts = 4;
        face[numFaces - 7].vert = new VertexID[face[numFaces - 7].nVerts];
        face[numFaces - 7].vert[0].vertIndex = numVerts - 18;
        face[numFaces - 7].vert[1].vertIndex = numVerts - 34;
        face[numFaces - 7].vert[2].vertIndex = numVerts - 33;
        face[numFaces - 7].vert[3].vertIndex = numVerts - 17;

        // right
        // front
        face[numFaces - 8].nVerts = 4;
        face[numFaces - 8].vert = new VertexID[face[numFaces - 8].nVerts];
        face[numFaces - 8].vert[0].vertIndex = numVerts - 32;
        face[numFaces - 8].vert[1].vertIndex = numVerts - 31;
        face[numFaces - 8].vert[2].vertIndex = numVerts - 11;
        face[numFaces - 8].vert[3].vertIndex = numVerts - 12;
        // back
        face[numFaces - 9].nVerts = 4;
        face[numFaces - 9].vert = new VertexID[face[numFaces - 9].nVerts];
        face[numFaces - 9].vert[0].vertIndex = numVerts - 20;
        face[numFaces - 9].vert[1].vertIndex = numVerts - 19;
        face[numFaces - 9].vert[2].vertIndex = numVerts - 35;
        face[numFaces - 9].vert[3].vertIndex = numVerts - 36;

        // circular faces
        for(int i=0; i<nSegments; i++) {
            // front
            face[i].nVerts = 4;
            face[i].vert = new VertexID[face[i].nVerts];
            face[i].vert[0].vertIndex = i % nSegments;
            face[i].vert[1].vertIndex = i % nSegments + 2*nSegments;
            face[i].vert[2].vertIndex = (i+1) % nSegments + 2*nSegments;
            face[i].vert[3].vertIndex = (i+1) % nSegments;
            // back
            face[i + nSegments].nVerts = 4;
            face[i + nSegments].vert = new VertexID[face[i + nSegments].nVerts];
            face[i + nSegments].vert[0].vertIndex = i % nSegments + nSegments;
            face[i + nSegments].vert[1].vertIndex = (i+1) % nSegments + nSegments;
            face[i + nSegments].vert[2].vertIndex = (i+1) % nSegments + 3*nSegments;
            face[i + nSegments].vert[3].vertIndex = i % nSegments + 3*nSegments;
        }
    }
};

void drawCoordinateSystem() {
    // x axis
    glColor3f(1,0,0);
    glBegin(GL_LINES);
    glVertex3f(0, 0, 0);
    glVertex3f(5, 0,0);
    glEnd();

    // y axis
    glColor3f(0,1,0);
    glBegin(GL_LINES);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 5, 0);
    glEnd();  

    // z axis
    glColor3f(0,0,1);
    glBegin(GL_LINES);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 5);
    glEnd(); 
}

// pillar
float pillarBodyWidth = 0.2;
float pillarBodyLength = 0.6;
float pillarBodyHeight = 3.8;
RectangularPrism pillarBody = RectangularPrism(pillarBodyWidth, pillarBodyLength, pillarBodyHeight);

float pillarBottomWidth = 1.5*pillarBodyLength;
float pillarBottomLength = pillarBodyHeight;
float pillarBottomHeight = pillarBodyWidth;
RectangularPrism pillarBottom = RectangularPrism(pillarBottomWidth, pillarBottomLength, pillarBottomHeight);

float pillarTopWidth = pillarBodyWidth;
float pillarTopLength = 1.1*pillarBottomLength;
float pillarTopHeight = 2*pillarBodyLength;
RectangularPrism pillarTop = RectangularPrism(pillarTopWidth, pillarTopLength, pillarTopHeight);

// rotation
float rotationRadius = 1.2*pillarBodyLength;
float rotationHeight = pillarBodyWidth;
float rotationSegments = 100;
Cylinder rotation = Cylinder(rotationRadius, rotationHeight, rotationSegments);

// crank
float crankWidth = 0.2;
float crankHeight = rotationHeight;
float crankLength = rotationRadius + crankWidth / 2;
float crankSegments = 50;
Crank crank = Crank(crankWidth, crankHeight, crankLength, crankSegments);

// bearing
float bearingBottomWidth = pillarBodyWidth;
float bearingBottomHeight= rotationHeight + crankHeight;
float bearingBottomLength = 0.8*pillarTopHeight;
float bearingBodyLength = 0.6*bearingBottomLength;
float bearingBodyHeight = 4*bearingBottomHeight;
float bearingRadius = bearingBodyLength / 3;
int bearingSegments = 100;
Bearing bearing = Bearing(bearingBottomWidth, bearingBottomHeight, bearingBottomLength,
    bearingBodyLength, bearingBodyLength, bearingRadius, bearingSegments);

// cylinder
float sliderRadius = bearingRadius;
float sliderHeight = 2*pillarTop.length;
float sliderSegments = 100;
Cylinder slider = Cylinder(sliderRadius, sliderHeight, sliderSegments);

// ubar
float ubarWidth = crankWidth;//pillarBodyWidth;
float ubarHeight = 0.35*pillarBodyHeight;
float ubarLength = 3*ubarWidth;
UBar ubar = UBar(ubarWidth, ubarHeight, ubarLength);

// nut
float nutSize = ubarWidth;
float nutRadius = 0.6*nutSize/2;
float nutHeight = ubarWidth;
int nutSegments = 100;
Nut nut = Nut(nutSize, nutHeight, nutRadius, nutSegments);

// linkbar
float linkbarLength = pillarBodyHeight;
float linkbarWidth = 3*crankWidth;
float linkbarHeight = crankHeight;
float linkbarSlotLength = 0.55*linkbarLength;
float linkBarSlotShift = 0.1*linkbarSlotLength;
float linkbarSlotWidth = crankWidth;
int linkbarSegments = 100;
LinkBar linkbar = LinkBar(linkbarLength, linkbarWidth, linkbarHeight, linkbarSlotLength, linkBarSlotShift, linkbarSlotWidth, linkbarSegments);

// bolt
float boltRadius = nutRadius;
float boltHeight = nutHeight;
float boltSegments = 100;
Cylinder bolt = Cylinder(boltRadius, boltHeight, boltSegments);

// first tile
float firstTileSize = 1.0;
int firstTileSegments = 25;
FirstTile firstTile = FirstTile(firstTileSize, firstTileSegments);

// second tile
float secondTileSize = firstTileSize;
SecondTile secondTile = SecondTile(secondTileSize);

// third tile
float thirdTileSize = 1.0;
int thirdTileSegments = 25;
ThirdTile thirdTile = ThirdTile(thirdTileSize, thirdTileSegments);


const int NUM_TILES_VERTICAL = 50;
const int NUM_TILES_HORIZONTAL = 50;

int** tiles;

void initTiles();
void initObjects() {
    pillarBody.init();
    pillarBottom.init();
    pillarTop.init();
    rotation.init();
    crank.init();
    bearing.init();
    slider.init();
    ubar.init();
    nut.init();
    linkbar.init();
    bolt.init();

    firstTile.init();
    secondTile.init();
    thirdTile.init();

    initTiles();
}

bool drawWireFrame = false;

// fixed
void drawPillarBody() {
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(-pillarBody.length / 2, pillarBottomHeight, -pillarBody.width / 2);


    if(!drawWireFrame) {
        pillarBody.setupMaterial(ColorArr[2]);
        pillarBody.Draw();
    } else {
        pillarBody.setColor(ColorArr[2]);
        pillarBody.DrawWireframe();
    }


    glPopMatrix();
}

// fixed
void drawPillarBottom() {
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(-pillarBottom.length / 2, 0, -pillarBottom.width / 2);

    if(!drawWireFrame) {
        pillarBottom.setupMaterial(ColorArr[2]);
        pillarBottom.Draw();
    } else {
        pillarBottom.setColor(ColorArr[2]);
        pillarBottom.DrawWireframe();
    }


    glPopMatrix();
}

// fixed
void drawPillarTop() {
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(-pillarTop.length / 2, pillarBottom.height + pillarBody.height, -pillarTop.width / 2);

    if(!drawWireFrame) {
        // pillarTop.DrawColor();
        pillarTop.setupMaterial(ColorArr[2]);
        pillarTop.Draw();
    } else {
        pillarTop.setColor(ColorArr[2]);
        pillarTop.DrawWireframe();
    }


    glPopMatrix();
}

int rotateAngle = 0; // 0: // -Ox -> rotate counter clockwise (Ox->Oy)
void drawRotation() {
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(0, pillarBottom.height + pillarBody.height / 2, rotation.height / 2 + pillarBody.width / 2);
    glRotatef(rotateAngle, 0, 0, 1);
    glRotatef(90, 1, 0, 0);
    glTranslatef(0, -rotation.height / 2, 0);


    if(!drawWireFrame) {
        rotation.setupMaterial(ColorArr[0]);
        rotation.Draw();
    } else {
        rotation.setColor(ColorArr[0]);
        rotation.DrawWireframe();
    }


    glPopMatrix();
} 

void drawCrank() {
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    // glScalef(3,3,3);
    glTranslatef(0, pillarBottom.height + pillarBody.height / 2, crank.height / 2 + pillarBody.width / 2 + rotation.height);
    glRotatef(rotateAngle, 0, 0, 1);
    glRotatef(90, 1, 0, 0);
    glTranslatef(-crank.length / 2 + crank.width / 2, -crank.height / 2, 0);
    // glScalef(3,3,3);


    if(!drawWireFrame) {
        crank.setupMaterial(ColorArr[9]);
        crank.Draw();
    } else {
        crank.setColor(ColorArr[9]);
        crank.DrawWireframe();
    }


    glPopMatrix();
}

// fixed
void drawBearing() {
    // left bearing
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(-0.9*pillarTop.length/2, pillarBottom.height + pillarBody.height + pillarTop.height/2, pillarTop.width/2);
    glRotatef(90, 1, 0, 0);
    glRotatef(90, 0, 1, 0);

    if(!drawWireFrame) {
        bearing.setupMaterial(ColorArr[0]);
        bearing.Draw();
    } else {
        bearing.setColor(ColorArr[0]);
        bearing.DrawWireframe();
    }


    glPopMatrix();


    // right bearing
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(0.9*pillarTop.length/2, pillarBottom.height + pillarBody.height + pillarTop.height/2, pillarTop.width/2);
    glRotatef(90, 1, 0, 0);
    glRotatef(90, 0, 1, 0);

    if(!drawWireFrame) {
        bearing.setupMaterial(ColorArr[0]);
        bearing.Draw();
    } else {
        bearing.setColor(ColorArr[0]);
        bearing.DrawWireframe();
    }

    glPopMatrix();
}

// fixed
void drawThirdBolt() {
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(0, pillarBottom.height + 0.1*pillarBody.height, pillarBody.width/2);
    glRotatef(90, 1, 0, 0);
    glScalef(2,1.2*(rotation.height + crank.height + linkbar.height)/bolt.height,2);


    if(!drawWireFrame) {
        bolt.setupMaterial(ColorArr[7]);
        bolt.Draw();
    } else {
        bolt.setColor(ColorArr[7]);
        bolt.DrawWireframe();
    }


    glPopMatrix();
}

float secondBoltPos = 0.5; // 0->1
void drawSecondBolt() {
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(0, pillarBottom.height + pillarBody.height/2, pillarBody.width/2 + rotation.height + crank.height);
    glRotatef(rotateAngle, 0, 0, 1);
    glTranslated(-secondBoltPos * (crank.length - crank.width), 0, 0); // 0 -> crank.length - crank.width
    glRotatef(90, 1, 0, 0);

    if(!drawWireFrame) {
        bolt.setupMaterial(ColorArr[7]);
        bolt.Draw();
    } else {
        bolt.setColor(ColorArr[7]);
        bolt.DrawWireframe();
    }


    glPopMatrix();
}

void drawSecondNut() {
    float d = 0.5*pillarBody.height - 0.1*pillarBody.height;
    float r = secondBoltPos * (crank.length - crank.width);
    float temp1 = r*sin((90 + rotateAngle) * PI / 180.0);
    float temp2 = r*cos((90 + rotateAngle) * PI / 180) + d;
    float thirdBoltAngle = atan(temp1 / temp2) * 180 / PI;


    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(0, pillarBottom.height + pillarBody.height/2, pillarBody.width/2 + rotation.height + crank.height);
    glRotatef(rotateAngle, 0, 0, 1); // rotate around rotation's center
    glTranslated(-r, 0, 0); // 0 -> crank.length - crank.width
    glRotatef(-rotateAngle, 0, 0, 1); // rotate around origin: to keep nut not rotate after glRotatef(rotateAngle, 0, 0, 1);
    glRotatef(thirdBoltAngle, 0, 0, 1); // rotate around origin: to make edges parallel with (secondBolt, thirdBolt) line


    if(!drawWireFrame) {
        nut.setupMaterial(ColorArr[0]);
        nut.Draw();
    } else {
        nut.setColor(ColorArr[0]);
        nut.DrawWireframe();
    }


    glPopMatrix();
}

void drawLinkBar() {
    float d = 0.5*pillarBody.height - 0.1*pillarBody.height;
    float r = secondBoltPos * (crank.length - crank.width);
    float temp1 = r*sin((90 + rotateAngle) * PI / 180.0);
    float temp2 = r*cos((90 + rotateAngle) * PI / 180) + d;
    float thirdBoltAngle = atan(temp1 / temp2) * 180 / PI;

    
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(0, pillarBottom.height + 0.1*pillarBody.height, 0.5*pillarBody.width + rotation.height + crank.height);
    glRotatef(thirdBoltAngle, 0, 0, 1);
    glRotatef(-90, 0, 0, 1);
    glTranslatef(-linkbar.length/2 + linkbar.width/2, 0, 0);


    if(!drawWireFrame) {
        linkbar.setupMaterial(ColorArr[1]);
        linkbar.Draw();
    } else {
        linkbar.setColor(ColorArr[1]);
        linkbar.DrawWireframe();
    }


    glPopMatrix();
}

void drawFirstBolt() {
    float d = 0.5*pillarBody.height - 0.1*pillarBody.height;
    float r = secondBoltPos * (crank.length - crank.width);
    float temp1 = r*sin((90 + rotateAngle) * PI / 180.0);
    float temp2 = r*cos((90 + rotateAngle) * PI / 180) + d;
    float thirdBoltAngle = atan(temp1 / temp2) * 180 / PI;

    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(0, pillarBottom.height + 0.1*pillarBody.height, pillarBody.width/2 + rotation.height + crank.height + linkbar.height);
    glRotatef(thirdBoltAngle, 0, 0, 1);
    glTranslatef(0, linkbar.length - linkbar.width, 0);
    glRotatef(90, 1, 0, 0);



    if(!drawWireFrame) {
        bolt.setupMaterial(ColorArr[7]);
        bolt.Draw();
    } else {
        bolt.setColor(ColorArr[7]);
        bolt.DrawWireframe();
    }


    glPopMatrix();
}

void drawFirstNut() {
    float d = 0.5*pillarBody.height - 0.1*pillarBody.height;
    float r = secondBoltPos * (crank.length - crank.width);
    float temp1 = r*sin((90 + rotateAngle) * PI / 180.0);
    float temp2 = r*cos((90 + rotateAngle) * PI / 180) + d;
    float thirdBoltAngle = atan(temp1 / temp2) * 180 / PI;

    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(0, pillarBottom.height + 0.1*pillarBody.height, pillarBody.width/2 + rotation.height + crank.height + linkbar.height);
    glRotatef(thirdBoltAngle, 0, 0, 1); // rotate around third bolt
    glTranslatef(0, linkbar.length - linkbar.width, 0);
    glRotatef(-thirdBoltAngle, 0, 0, 1); // rotate around origin -> no rotate after glRotatef(thirdBoltAngle, 0, 0, 1);


    if(!drawWireFrame) {
        nut.setupMaterial(ColorArr[0]);
        nut.Draw();
    } else {
        nut.setColor(ColorArr[0]);
        nut.DrawWireframe();
    }


    glPopMatrix();
}

void drawSlider() {
    float d = 0.5*pillarBody.height - 0.1*pillarBody.height;
    float r = secondBoltPos * (crank.length - crank.width);
    float temp1 = r*sin((90 + rotateAngle) * PI / 180.0);
    float temp2 = r*cos((90 + rotateAngle) * PI / 180) + d;
    float thirdBoltAngle = atan(temp1 / temp2) * 180 / PI;

    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef((linkbar.length - linkbar.width)*sin(-thirdBoltAngle * PI / 180), pillarBottom.height + pillarBody.height + pillarTop.height/2, pillarTop.width/2 + bearing.bottomHeight + bearing.bodyHeight - bearing.bodyLength/2);
    glRotatef(90, 0, 0, 1);
    glTranslatef(0, -slider.height/2, 0);

    if(!drawWireFrame) {
        slider.setupMaterial(ColorArr[6]);
        slider.Draw();
    } else {
        slider.setColor(ColorArr[6]);
        slider.DrawWireframe();
    }


    glPopMatrix();
}

void drawUBar() {
    float d = 0.5*pillarBody.height - 0.1*pillarBody.height;
    float r = secondBoltPos * (crank.length - crank.width);
    float temp1 = r*sin((90 + rotateAngle) * PI / 180.0);
    float temp2 = r*cos((90 + rotateAngle) * PI / 180) + d;
    float thirdBoltAngle = atan(temp1 / temp2) * 180 / PI;

    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef((linkbar.length - linkbar.width)*sin(-thirdBoltAngle * PI / 180), pillarBottom.height + pillarBody.height + pillarTop.height/2 - slider.radius, pillarTop.width/2 + bearing.bottomHeight + bearing.bodyHeight - bearing.bodyLength/2);
    glRotatef(180, 1, 0, 0);
    glTranslatef(0, 0, -ubar.width/2);


    if(!drawWireFrame) {
        ubar.setupMaterial(ColorArr[6]);
        ubar.Draw();
    } else {
        ubar.setColor(ColorArr[6]);
        ubar.DrawWireframe();
    }


    glPopMatrix();
}



// draw floor
void initTiles() {
    // initialize random seed
    srand (time(NULL));

    // inittialize ramdom tiles
    tiles = new int*[NUM_TILES_HORIZONTAL];
    for(int i=0; i<NUM_TILES_HORIZONTAL; i++) {
        tiles[i] = new int[NUM_TILES_VERTICAL];
    }

    for(int i=0; i<NUM_TILES_HORIZONTAL; i++) {
        for(int j=0; j<NUM_TILES_VERTICAL; j++) {
            tiles[i][j] = rand() % 4 + 1; // generate random number from 1 to 4
        }
    }
}

void drawFirstTile(float x, float y) {
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(x, 0, y);
    glRotatef(90, 1, 0, 0);

    firstTile.setupFirstColor(ColorArr[11]);
    firstTile.setupSecondColor(ColorArr[7]);
    firstTile.DrawColor();
    
    glPopMatrix();
}

void drawSecondTile(float x, float y) {
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(x, 0, y);
    glRotatef(90, 1, 0, 0);

    secondTile.setupFirstColor(ColorArr[11]);
    secondTile.setupSecondColor(ColorArr[7]);
    // secondTile.DrawWireframe();
    secondTile.DrawColor();

    glPopMatrix();
}

void drawThirdTile(float x, float y) {
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(x, 0, y);
    glRotatef(90, 1, 0, 0);

    thirdTile.setupFirstColor(ColorArr[11]);
    thirdTile.setupSecondColor(ColorArr[7]);
    thirdTile.DrawColor();

    glPopMatrix();
}

void drawFourthTile(float x, float y) {
    glPushMatrix();

    glMatrixMode(GL_MODELVIEW);
    glTranslatef(x, 0, y);
    glRotatef(90, 1, 0, 0);
    glRotatef(-90, 0, 0, 1);

    firstTile.setupFirstColor(ColorArr[11]);
    firstTile.setupSecondColor(ColorArr[7]);
    firstTile.DrawColor();
    
    glPopMatrix();
}

void drawTile(int type, float x, float y) {
    switch (type)
    {
    case 1:
        drawFirstTile(x,y);
        break;
    case 2:
        drawSecondTile(x,y);
        break;
    case 3:
        drawThirdTile(x,y);
        break;
    case 4:
        drawFourthTile(x,y);
        break;    
    default:
        break;
    }
}

void drawFloor() {
    for(int i=-NUM_TILES_HORIZONTAL/2; i<NUM_TILES_HORIZONTAL/2; i++) {
        for(int j=-NUM_TILES_VERTICAL/2; j<NUM_TILES_VERTICAL/2; j++) {
            float x = firstTileSize * i;
            float y = firstTileSize * j;
            drawTile(tiles[i+NUM_TILES_HORIZONTAL/2][j+NUM_TILES_VERTICAL/2], x, y);
        }
    }
}


void drawObjects() {
    drawPillarBody();
    drawPillarBottom();
    drawPillarTop();
    drawRotation();
    drawCrank();
    drawBearing();
    drawSlider();
    drawUBar();
    drawThirdBolt();
    drawSecondBolt();
    drawSecondNut();
    drawLinkBar();
    drawFirstBolt();
    drawFirstNut();

    drawFloor();
}

float cameraAngle = -25;
float cameraHeight = 4;
float cameraDistance = 10;

bool secondLightOn = false;
bool perspective = true;

void display() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if(perspective) {
        // gluPerspective(75, 1, 1, 50);
        glFrustum(-1, 1, -0.5*SCREEN_HEIGHT/SCREEN_WIDTH, 1.5*SCREEN_HEIGHT/SCREEN_WIDTH, 1, 50.0);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(cameraDistance*sin(cameraAngle * PI / 180), cameraHeight, cameraDistance*cos(cameraAngle * PI / 180), 0, 0, 0, 0, pillarBottom.height + pillarBody.height/2, 0);
    } else {
        glOrtho(-10, 10, -10*SCREEN_HEIGHT/SCREEN_WIDTH, 10*SCREEN_HEIGHT/SCREEN_WIDTH, -1000, 1000);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        // cameraDistance = 0;
        gluLookAt(0, 20, 0, 0, 0, 0, 0, 0, 1);
    }



    //    gluLookAt(0, 0, 6, 0, 0, 0, 0, 1, 0);

    if(secondLightOn) {
        // set up first light source
        glEnable(GL_LIGHT1); 
        GLfloat	lightDiffuse[] = {0.8f, 0.8f, 0.8f, 1.0f};
        GLfloat	lightSpecular[] = {0.8f, 0.8f, 0.8f, 1.0f};
        GLfloat	lightAmbient[] = {0.5f, 0.5f, 0.5f, 1.0f};
        GLfloat light_position2[] = {-6.0f, 0.0f, 0.5f, 0.0f};
        glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
        glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiffuse);
        glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmbient);
        glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpecular);
    } else {
        glDisable(GL_LIGHT1);
    }

    if(drawWireFrame) {
        glDisable(GL_LIGHTING);
    } else {
        glEnable(GL_LIGHTING);
    }


    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
    
    // drawCoordinateSystem();

    drawObjects();
        
    glFlush();
    glutSwapBuffers();
}


void keyPressed(unsigned char key, int x, int y) {
    switch (key)
    {
    case '1':
        rotateAngle = (rotateAngle + 10) % 360; 
        break;
    case '2':
        rotateAngle = (rotateAngle - 10) % 360; 
        break;
    case '3':
        secondBoltPos = secondBoltPos - 0.1;
        if(secondBoltPos < 0) secondBoltPos = 0; 
        break;
    case '4':
        secondBoltPos = secondBoltPos + 0.1;
        if(secondBoltPos > 1) secondBoltPos = 1; 
        break;
    case '+':
        cameraDistance -= 0.1;
        break;
    case '-':
        cameraDistance += 0.1;
        break;
    case 'w':
        drawWireFrame = !drawWireFrame;
        break;
    case 'W':
        drawWireFrame = !drawWireFrame;
        break;
    case 'd':
        secondLightOn = !secondLightOn;
        break;
    case 'D':
        secondLightOn = !secondLightOn;
        break;
    case 'v':
        perspective = !perspective;
        break;
    case 'V':
        perspective = !perspective;
        break;
    default:
        break;
    }
    glutPostRedisplay();
}

void specialKeyPressed(int key, int x, int y) {
    switch (key)
    {
    case GLUT_KEY_UP:
        cameraHeight += 0.1;
        break;
    case GLUT_KEY_DOWN:
        cameraHeight -= 0.1;
        break;    
    case GLUT_KEY_LEFT:
        cameraAngle += 10;
        break;
    case GLUT_KEY_RIGHT:
        cameraAngle -= 10;
        break;  
    default:
        break;
    }
    glutPostRedisplay();
}


void init()
{
    float size = 10;

    glClearColor(1, 1, 1, 1);

    glFrontFace(GL_CCW);
    glEnable(GL_DEPTH_TEST); // enable depth calculation
    glEnable(GL_LIGHTING); // enable light mode

    // set up first light source
    glEnable(GL_LIGHT0); 
    GLfloat	lightDiffuse[] = {0.8f, 0.8f, 0.8f, 1.0f};
	GLfloat	lightSpecular[] = {0.8f, 0.8f, 0.8f, 1.0f};
	GLfloat	lightAmbient[] = {0.5f, 0.5f, 0.5f, 1.0f};
	GLfloat light_position1[] = {6.0f, 3.0f, 6.0f, 0.0f};
	glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);


    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // glOrtho(-0.5*size, 0.5*size, -0.5*size, 0.5*size, -1000, 1000);
    // glFrustum(-5, 5, -2.0*SCREEN_HEIGHT/SCREEN_WIDTH, 7.*SCREEN_HEIGHT/SCREEN_WIDTH, 1, 50.0);

    // gluPerspective(75, 1, 1, 50);
    
    
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void printInstruction() {
    cout << "1\t: Xoay ban xoay nguoc chieu kim dong ho\n";
    cout << "2\t: Xoay ban xoay cung chieu kim dong ho\n";
    cout << "3, 4\t: Dieu chinh vi tri cua chot 2\n";
    cout << "W, w\t: Chuyen doi qua lai giua che do khung day la to mau\n";
    cout << "V, v\t: Chuyen doi qua lai giua hai che do nhin khac nhau\n";
    cout << "D, d\t: Bat/tat nguon sang thu hai\n";
    cout << "+\t: Tang khoang cach camera\n";
    cout << "-\t: Giam khoang cach camera\n";
    cout << "up arrow\t: Tang chieu cao camera\n";
    cout << "down arrow\t: Giam chieu cao camera\n";
    cout << "<-\t: Quay camera theo chieu kim dong ho\n";
    cout << "->\t: Quay camera nguoc chieu kim dong ho\n";
}

int main(int argc, char** argv) {
    glutInit(&argc, argv); //initialize the tool kit
    glutInitDisplayMode(GLUT_DOUBLE |GLUT_RGB | GLUT_DEPTH);//set the display mode
    glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT); //set window size
    glutInitWindowPosition(100, 100); // set window position on screen
    glutCreateWindow("Bui Ngo Hoang Long - 1810283"); // open the screen window

    printInstruction();

    init();

    initObjects();

    glutDisplayFunc(display);
    glutKeyboardFunc(keyPressed);
    glutSpecialFunc(specialKeyPressed);
      
    glutMainLoop();

    for(int i=0; i<NUM_TILES_HORIZONTAL; i++) {
        delete[] tiles[i];
    }
    delete[] tiles;

    return 0;

}