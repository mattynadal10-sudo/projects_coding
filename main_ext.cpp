#define SOQT_NOT_DLL

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/SbTime.h>
#include <Inventor/SbVec2f.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/nodes/SoComplexity.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoNurbsCurve.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoTexture2.h>
#include <Inventor/nodes/SoTextureCoordinate2.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/sensors/SoTimerSensor.h>

#include <QApplication>
#include <QEvent>
#include <QKeyEvent>
#include <QObject>
#include <QWidget>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>

static const int N = 13;

// Units
static const float CM_PER_M = 100.0f;

// Physical parameters
static const float MASS = 2.2f;            // kg
static const float K = 18.0f;              // N/cm
static const float B = 2.4f;               // N*s/cm
static const float AIR_DRAG = 0.30f;       // N*s/cm
static const float G = 9.81f;              // m/s^2
static const float TOTAL_LENGTH = 115.0f;  // cm
static const float REST_LENGTH = TOTAL_LENGTH / (N - 1);

// Simulation parameters
static const float DT = 0.0025f;
static const int SUBSTEPS = 2;
static const float PARTICLE_RADIUS = 2.0f;

// Constraint between q5 and q6
static const int C_A = 4;
static const int C_B = 5;
static float gConstraintLength = 0.0f;

static const float BAUMGARTE_BETA = 3.5f;

// Wall
static const float WALL_X = 42.0f;
static const float WALL_RESTITUTION = 0.82f;
static const float WALL_TANGENTIAL_DAMP = 0.999f;
static const float WALL_EPS = 0.6f;
static const float WALL_SEPARATING_SPEED = 18.0f; // cm/s

struct Particle
{
    SbVec3f pos;   // cm
    SbVec3f vel;   // cm/s
    bool fixed;
};

static std::vector<Particle> gParticles(N);
static bool gWindOn = false;
static float gSimTime = 0.0f;
static float gLastLambda = 0.0f;

// Scene refs
static SoCoordinate3* gCurveCoords = nullptr;
static std::vector<SoTranslation*> gParticleTranslations;
static SoTimerSensor* gTimer = nullptr;

// Logging
static std::ofstream gDistanceLog;

float constrainedDistance()
{
    return (gParticles[C_A].pos - gParticles[C_B].pos).length();
}

void logDistance()
{
    if (gDistanceLog.is_open())
    {
        gDistanceLog << constrainedDistance() << "\n";
        gDistanceLog.flush();
    }
}

void initParticles()
{
    const float anchorX = -8.0f;
    const float anchorY = 155.0f;
    const float anchorZ = 0.0f;

    gParticles[0].pos = SbVec3f(anchorX, anchorY, anchorZ);
    gParticles[0].vel = SbVec3f(0.0f, 0.0f, 0.0f);
    gParticles[0].fixed = true;

    float currentY = anchorY;

    for (int i = 1; i < N; ++i)
    {
        int massesBelow = N - i;
        float supportedWeight = massesBelow * MASS * G; // N
        float extension = supportedWeight / K;          // cm
        float segmentLength = REST_LENGTH + extension;  // cm

        currentY -= segmentLength;

        gParticles[i].pos = SbVec3f(anchorX, currentY, anchorZ);
        gParticles[i].vel = SbVec3f(0.0f, 0.0f, 0.0f);
        gParticles[i].fixed = false;
    }

    gConstraintLength = (gParticles[C_A].pos - gParticles[C_B].pos).length();

    gSimTime = 0.0f;
    gLastLambda = 0.0f;
    gWindOn = false;
}

SbVec3f windForce(int i)
{
    if (!gWindOn)
        return SbVec3f(0.0f, 0.0f, 0.0f);

    float fx = 65.0f + 5.0f * static_cast<float>(i);
    return SbVec3f(fx, 0.0f, 0.0f);
}

void applySpringDamper(int i, int j, std::vector<SbVec3f>& forces)
{
    SbVec3f d = gParticles[j].pos - gParticles[i].pos;
    float L = d.length();
    if (L < 1e-6f)
        return;

    SbVec3f dir = d / L;

    float springMag = K * (L - REST_LENGTH);

    SbVec3f relVel = gParticles[j].vel - gParticles[i].vel;
    float relAlong = relVel.dot(dir);
    float damperMag = B * relAlong;

    SbVec3f f = (springMag + damperMag) * dir;

    if (!gParticles[i].fixed) forces[i] += f;
    if (!gParticles[j].fixed) forces[j] -= f;
}

float solveLambda(const std::vector<SbVec3f>& forcesForPrediction)
{
    const Particle& pa = gParticles[C_A];
    const Particle& pb = gParticles[C_B];

    SbVec3f d = pa.pos - pb.pos;
    float d2 = d.dot(d);

    if (d2 < 1e-10f)
        return 0.0f;

    SbVec3f vaStar = pa.vel;
    SbVec3f vbStar = pb.vel;

    if (!pa.fixed) vaStar += DT * (CM_PER_M / MASS) * forcesForPrediction[C_A];
    if (!pb.fixed) vbStar += DT * (CM_PER_M / MASS) * forcesForPrediction[C_B];

    SbVec3f dvStar = vaStar - vbStar;

    float JvStar = 2.0f * d.dot(dvStar);
    float Cval = d2 - gConstraintLength * gConstraintLength;

    float invMa = pa.fixed ? 0.0f : (CM_PER_M / MASS);
    float invMb = pb.fixed ? 0.0f : (CM_PER_M / MASS);

    float JMInvJT = 4.0f * d2 * (invMa + invMb);
    if (JMInvJT < 1e-12f)
        return 0.0f;

    float rhs = JvStar + (BAUMGARTE_BETA / DT) * Cval;
    return -rhs / (DT * JMInvJT);
}

void applyConstraintToVelocities(float lambda)
{
    SbVec3f d = gParticles[C_A].pos - gParticles[C_B].pos;

    SbVec3f Ja =  2.0f * d;
    SbVec3f Jb = -2.0f * d;

    if (!gParticles[C_A].fixed)
        gParticles[C_A].vel += DT * (CM_PER_M / MASS) * (Ja * lambda);

    if (!gParticles[C_B].fixed)
        gParticles[C_B].vel += DT * (CM_PER_M / MASS) * (Jb * lambda);
}

void projectConstraintPosition()
{
    SbVec3f d = gParticles[C_A].pos - gParticles[C_B].pos;
    float L = d.length();

    if (L < 1e-8f)
        return;

    float error = L - gConstraintLength;
    if (std::fabs(error) < 1e-6f)
        return;

    SbVec3f dir = d / L;

    float wa = gParticles[C_A].fixed ? 0.0f : 1.0f / MASS;
    float wb = gParticles[C_B].fixed ? 0.0f : 1.0f / MASS;
    float wsum = wa + wb;

    if (wsum < 1e-12f)
        return;

    if (!gParticles[C_A].fixed)
        gParticles[C_A].pos -= dir * (error * wa / wsum);

    if (!gParticles[C_B].fixed)
        gParticles[C_B].pos += dir * (error * wb / wsum);
}

void resolveWallCollisionForParticle(int i)
{
    if (gParticles[i].fixed)
        return;

    const float wallLimit = WALL_X - PARTICLE_RADIUS - WALL_EPS;

    if (gParticles[i].pos[0] > wallLimit)
    {
        float penetration = gParticles[i].pos[0] - wallLimit;

        gParticles[i].pos[0] = wallLimit - penetration;

        if (gParticles[i].vel[0] > 0.0f)
            gParticles[i].vel[0] = -WALL_RESTITUTION * gParticles[i].vel[0];

        if (gParticles[i].vel[0] > -WALL_SEPARATING_SPEED)
            gParticles[i].vel[0] = -WALL_SEPARATING_SPEED;

        gParticles[i].vel[1] *= WALL_TANGENTIAL_DAMP;
        gParticles[i].vel[2] *= WALL_TANGENTIAL_DAMP;
    }
}

void resolveWallCollisions()
{
    for (int i = 0; i < N; ++i)
        resolveWallCollisionForParticle(i);

    projectConstraintPosition();
}

void simulateStep()
{
    std::vector<SbVec3f> forces(N, SbVec3f(0.0f, 0.0f, 0.0f));

    for (int i = 0; i < N; ++i)
    {
        if (gParticles[i].fixed)
            continue;

        forces[i] += SbVec3f(0.0f, -MASS * G, 0.0f);
        forces[i] += -AIR_DRAG * gParticles[i].vel;
        forces[i] += windForce(i);
    }

    for (int i = 0; i < N - 1; ++i)
        applySpringDamper(i, i + 1, forces);

    for (int i = 0; i < N; ++i)
    {
        if (gParticles[i].fixed)
            continue;

        SbVec3f acc = (CM_PER_M / MASS) * forces[i];
        gParticles[i].vel += DT * acc;
    }

    gLastLambda = solveLambda(forces);
    applyConstraintToVelocities(gLastLambda);

    for (int i = 0; i < N; ++i)
    {
        if (gParticles[i].fixed)
            continue;

        gParticles[i].pos += DT * gParticles[i].vel;
    }

    projectConstraintPosition();
    resolveWallCollisions();

    gSimTime += DT;
    logDistance();
}

void updateScene()
{
    for (int i = 0; i < N; ++i)
        gParticleTranslations[i]->translation.setValue(gParticles[i].pos);

    const int numCtrl = N + 4;
    SbVec3f ctrl[numCtrl];

    ctrl[0] = gParticles[0].pos;
    ctrl[1] = gParticles[0].pos;

    for (int i = 0; i < N; ++i)
        ctrl[i + 2] = gParticles[i].pos;

    ctrl[numCtrl - 2] = gParticles[N - 1].pos;
    ctrl[numCtrl - 1] = gParticles[N - 1].pos;

    gCurveCoords->point.setValues(0, numCtrl, ctrl);
}

void timerCB(void*, SoSensor*)
{
    for (int i = 0; i < SUBSTEPS; ++i)
        simulateStep();

    updateScene();
}

SoSeparator* makeParticleNode(SoTranslation*& trOut)
{
    SoSeparator* sep = new SoSeparator;

    SoMaterial* mat = new SoMaterial;
    mat->diffuseColor.setValue(0.88f, 0.18f, 0.18f);
    sep->addChild(mat);

    trOut = new SoTranslation;
    sep->addChild(trOut);

    SoSphere* s = new SoSphere;
    s->radius = PARTICLE_RADIUS;
    sep->addChild(s);

    return sep;
}

SoSeparator* makeWallNode()
{
    SoSeparator* wallSep = new SoSeparator;

    SoTexture2* tex = new SoTexture2;
    tex->filename = "wall.jpg";
    wallSep->addChild(tex);

    SoMaterial* mat = new SoMaterial;
    mat->diffuseColor.setValue(1.0f, 1.0f, 1.0f);
    wallSep->addChild(mat);

    SoCoordinate3* coords = new SoCoordinate3;
    coords->point.set1Value(0, SbVec3f(WALL_X, -80.0f, -110.0f));
    coords->point.set1Value(1, SbVec3f(WALL_X, -80.0f,  110.0f));
    coords->point.set1Value(2, SbVec3f(WALL_X,  220.0f, 110.0f));
    coords->point.set1Value(3, SbVec3f(WALL_X,  220.0f,-110.0f));
    wallSep->addChild(coords);

    SoTextureCoordinate2* tc = new SoTextureCoordinate2;
    tc->point.set1Value(0, SbVec2f(0.0f, 0.0f));
    tc->point.set1Value(1, SbVec2f(1.0f, 0.0f));
    tc->point.set1Value(2, SbVec2f(1.0f, 1.0f));
    tc->point.set1Value(3, SbVec2f(0.0f, 1.0f));
    wallSep->addChild(tc);

    SoFaceSet* face = new SoFaceSet;
    face->numVertices.set1Value(0, 4);
    wallSep->addChild(face);

    return wallSep;
}

class KeyFilter : public QObject
{
public:
    explicit KeyFilter(QObject* parent = nullptr) : QObject(parent) {}

protected:
    bool eventFilter(QObject*, QEvent* event) override
    {
        if (event->type() == QEvent::KeyPress)
        {
            QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);

            if (keyEvent->key() == Qt::Key_W)
            {
                gWindOn = !gWindOn;
                std::printf("Wind %s\n", gWindOn ? "ON" : "OFF");
                std::fflush(stdout);
                return true;
            }

            if (keyEvent->key() == Qt::Key_R)
            {
                initParticles();
                updateScene();
                std::printf("Simulation reset\n");
                std::fflush(stdout);
                return true;
            }

            if (keyEvent->key() == Qt::Key_Escape)
            {
                QApplication::quit();
                return true;
            }
        }

        return false;
    }
};

int main(int argc, char** argv)
{
    QWidget* window = SoQt::init(argv[0]);
    if (!window) return 1;

    gDistanceLog.open("distance.txt", std::ios::out | std::ios::app);

    initParticles();

    SoSeparator* root = new SoSeparator;
    root->ref();

    SoPerspectiveCamera* cam = new SoPerspectiveCamera;
    cam->position.setValue(110.0f, 65.0f, 280.0f);
    cam->nearDistance = 1.0f;
    cam->farDistance = 2500.0f;
    root->addChild(cam);

    SoDirectionalLight* light = new SoDirectionalLight;
    light->direction.setValue(-0.45f, -0.55f, -1.0f);
    root->addChild(light);

    root->addChild(makeWallNode());

    SoSeparator* ropeSep = new SoSeparator;

    SoMaterial* ropeMat = new SoMaterial;
    ropeMat->diffuseColor.setValue(0.10f, 0.24f, 0.88f);
    ropeSep->addChild(ropeMat);

    SoDrawStyle* ropeStyle = new SoDrawStyle;
    ropeStyle->style = SoDrawStyle::LINES;
    ropeStyle->lineWidth = 5.0f;
    ropeSep->addChild(ropeStyle);

    SoComplexity* comp = new SoComplexity;
    comp->value = 0.85f;
    ropeSep->addChild(comp);

    gCurveCoords = new SoCoordinate3;
    ropeSep->addChild(gCurveCoords);

    SoNurbsCurve* curve = new SoNurbsCurve;
    const int numCtrl = N + 4;
    curve->numControlPoints = numCtrl;

    float knots[numCtrl + 4];
    for (int i = 0; i < 4; ++i) knots[i] = 0.0f;
    for (int i = 4; i < numCtrl; ++i) knots[i] = static_cast<float>(i - 3);
    float lastKnot = static_cast<float>(numCtrl - 3);
    for (int i = numCtrl; i < numCtrl + 4; ++i) knots[i] = lastKnot;

    curve->knotVector.setValues(0, numCtrl + 4, knots);
    ropeSep->addChild(curve);

    root->addChild(ropeSep);

    gParticleTranslations.resize(N);
    for (int i = 0; i < N; ++i)
    {
        SoTranslation* tr = nullptr;
        SoSeparator* p = makeParticleNode(tr);
        gParticleTranslations[i] = tr;
        root->addChild(p);
    }

    SoQtExaminerViewer* viewer = new SoQtExaminerViewer(window);
    viewer->setSceneGraph(root);
    viewer->setTitle("Constrained Rope With Strong Wall Crash");
    viewer->setBackgroundColor(SbColor(0.93f, 0.93f, 0.97f));
    viewer->show();

    updateScene();
    cam->viewAll(root, viewer->getViewportRegion());

    KeyFilter* filter = new KeyFilter(window);
    window->installEventFilter(filter);
    window->setFocusPolicy(Qt::StrongFocus);
    window->setFocus();

    std::printf("\nControls:\n");
    std::printf("  W   -> toggle wind ON/OFF\n");
    std::printf("  R   -> reset rope\n");
    std::printf("  Esc -> quit\n");
    std::printf("  Logging to distance.txt\n\n");
    std::fflush(stdout);

    gTimer = new SoTimerSensor(timerCB, nullptr);
    gTimer->setInterval(SbTime(1.0 / 60.0));
    gTimer->schedule();

    SoQt::show(window);
    SoQt::mainLoop();

    if (gTimer)
    {
        gTimer->unschedule();
        delete gTimer;
        gTimer = nullptr;
    }

    if (gDistanceLog.is_open())
        gDistanceLog.close();

    root->unref();
    return 0;
}