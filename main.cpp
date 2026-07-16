#define SOQT_NOT_DLL

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>

#include <Inventor/SbTime.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/events/SoKeyboardEvent.h>
#include <Inventor/nodes/SoComplexity.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoNurbsCurve.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/sensors/SoTimerSensor.h>

#include <QApplication>
#include <QWidget>
#include <Qt>

#include <cmath>
#include <cstdio>
#include <vector>

static const int N = 13;

// Assignment data
static const float MASS = 1.5f;          // kg
static const float K = 4.6f;             // N/cm
static const float B = 2.2f;             // N*s/cm
static const float AIR_DRAG = 0.9f;      // linear drag coefficient
static const float G = 9.81f;            // m/s^2
static const float TOTAL_LENGTH = 30.0f; // cm
static const float REST_LENGTH = TOTAL_LENGTH / (N - 1);

// Simulation parameters
static const float DT = 0.004f;          // s
static const int SUBSTEPS = 2;
static const float PARTICLE_RADIUS = 0.7f;

struct Particle
{
    SbVec3f pos;   // cm
    SbVec3f vel;   // cm/s
    bool fixed;
};

static std::vector<Particle> gParticles(N);
static bool gWindOn = false;
static float gTime = 0.0f;

// Scene references
static SoCoordinate3* gCurveCoords = nullptr;
static std::vector<SoTranslation*> gParticleTranslations;
static SoTimerSensor* gTimer = nullptr;

// ------------------------------
// Initialize rope in static equilibrium
// ------------------------------
void initParticles()
{
    gParticles[0].pos = SbVec3f(0.0f, 50.0f, 0.0f);
    gParticles[0].vel = SbVec3f(0.0f, 0.0f, 0.0f);
    gParticles[0].fixed = true;

    float currentY = 50.0f;

    for (int i = 1; i < N; ++i)
    {
        // Spring (i-1) supports particles i ... N-1
        int massesBelow = N - i;
        float supportedWeight = massesBelow * MASS * G; // N
        float extension = supportedWeight / K;          // cm
        float segmentLength = REST_LENGTH + extension;  // cm

        currentY -= segmentLength;

        gParticles[i].pos = SbVec3f(0.0f, currentY, 0.0f);
        gParticles[i].vel = SbVec3f(0.0f, 0.0f, 0.0f);
        gParticles[i].fixed = false;
    }

    gTime = 0.0f;
}

// ------------------------------
// Wind force in +x and -z
// ------------------------------
SbVec3f windForce(int i)
{
    if (!gWindOn)
        return SbVec3f(0.0f, 0.0f, 0.0f);

    float phase = 0.20f * static_cast<float>(i);
    float fx = 40.0f + 18.0f * std::sin(2.2f * gTime + phase);
    float fz = -28.0f - 12.0f * std::cos(1.7f * gTime + 0.4f * phase);

    return SbVec3f(fx, 0.0f, fz);
}

// ------------------------------
// Spring-damper force
// ------------------------------
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

    if (!gParticles[i].fixed)
        forces[i] += f;
    if (!gParticles[j].fixed)
        forces[j] -= f;
}

// ------------------------------
// Simulation step
// ------------------------------
void simulateStep()
{
    std::vector<SbVec3f> forces(N, SbVec3f(0.0f, 0.0f, 0.0f));

    for (int i = 0; i < N; ++i)
    {
        if (gParticles[i].fixed)
            continue;

        forces[i] += SbVec3f(0.0f, -MASS * G, 0.0f);   // gravity
        forces[i] += -AIR_DRAG * gParticles[i].vel;    // drag
        forces[i] += windForce(i);                     // wind
    }

    for (int i = 0; i < N - 1; ++i)
        applySpringDamper(i, i + 1, forces);

    for (int i = 0; i < N; ++i)
    {
        if (gParticles[i].fixed)
            continue;

        SbVec3f acc = (1.0f / MASS) * forces[i];
        gParticles[i].vel += DT * acc;
        gParticles[i].pos += DT * gParticles[i].vel;
    }

    gTime += DT;
}

// ------------------------------
// Update scene
// ------------------------------
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

// ------------------------------
// Timer callback
// ------------------------------
void timerCB(void*, SoSensor*)
{
    for (int i = 0; i < SUBSTEPS; ++i)
        simulateStep();

    updateScene();
}

// ------------------------------
// Keyboard callback
// ------------------------------
void keyboardCB(void*, SoEventCallback* eventCB)
{
    const SoEvent* ev = eventCB->getEvent();

    if (!SoKeyboardEvent::isKeyPressEvent(ev, SoKeyboardEvent::ANY))
        return;

    const SoKeyboardEvent* kev = static_cast<const SoKeyboardEvent*>(ev);
    SoKeyboardEvent::Key key = kev->getKey();

    if (key == SoKeyboardEvent::W)
    {
        gWindOn = !gWindOn;
        std::printf("Wind %s\n", gWindOn ? "ON" : "OFF");
        std::fflush(stdout);
        eventCB->setHandled();
    }
    else if (key == SoKeyboardEvent::R)
    {
        initParticles();
        updateScene();
        std::printf("Simulation reset\n");
        std::fflush(stdout);
        eventCB->setHandled();
    }
    else if (key == SoKeyboardEvent::ESCAPE)
    {
        std::printf("Exiting...\n");
        std::fflush(stdout);
        QApplication::quit();
        eventCB->setHandled();
    }
}

// ------------------------------
// Create one particle node
// ------------------------------
SoSeparator* makeParticleNode(SoTranslation*& trOut)
{
    SoSeparator* sep = new SoSeparator;

    SoMaterial* mat = new SoMaterial;
    mat->diffuseColor.setValue(0.85f, 0.15f, 0.15f);
    sep->addChild(mat);

    trOut = new SoTranslation;
    sep->addChild(trOut);

    SoSphere* s = new SoSphere;
    s->radius = PARTICLE_RADIUS;
    sep->addChild(s);

    return sep;
}

int main(int argc, char** argv)
{
    QWidget* window = SoQt::init(argv[0]);
    if (!window) return 1;

    initParticles();

    SoSeparator* root = new SoSeparator;
    root->ref();

    // Camera
    SoPerspectiveCamera* cam = new SoPerspectiveCamera;
    cam->position.setValue(0.0f, 10.0f, 140.0f);
    cam->nearDistance = 1.0f;
    cam->farDistance = 1000.0f;
    root->addChild(cam);

    // Light
    SoDirectionalLight* light = new SoDirectionalLight;
    light->direction.setValue(-0.2f, -0.8f, -1.0f);
    root->addChild(light);

    // Rope curve
    SoSeparator* ropeSep = new SoSeparator;

    SoMaterial* ropeMat = new SoMaterial;
    ropeMat->diffuseColor.setValue(0.1f, 0.25f, 0.85f);
    ropeSep->addChild(ropeMat);

    SoDrawStyle* ropeStyle = new SoDrawStyle;
    ropeStyle->style = SoDrawStyle::LINES;
    ropeStyle->lineWidth = 3.0f;
    ropeSep->addChild(ropeStyle);

    SoComplexity* comp = new SoComplexity;
    comp->value = 0.8f;
    ropeSep->addChild(comp);

    gCurveCoords = new SoCoordinate3;
    ropeSep->addChild(gCurveCoords);

    SoNurbsCurve* curve = new SoNurbsCurve;
    const int numCtrl = N + 4;
    curve->numControlPoints = numCtrl;

    float knots[numCtrl + 4];
    for (int i = 0; i < 4; ++i)
        knots[i] = 0.0f;
    for (int i = 4; i < numCtrl; ++i)
        knots[i] = static_cast<float>(i - 3);
    float lastKnot = static_cast<float>(numCtrl - 3);
    for (int i = numCtrl; i < numCtrl + 4; ++i)
        knots[i] = lastKnot;

    curve->knotVector.setValues(0, numCtrl + 4, knots);
    ropeSep->addChild(curve);

    root->addChild(ropeSep);

    // Particle spheres
    gParticleTranslations.resize(N);
    for (int i = 0; i < N; ++i)
    {
        SoTranslation* tr = nullptr;
        SoSeparator* p = makeParticleNode(tr);
        gParticleTranslations[i] = tr;
        root->addChild(p);
    }

    // Keyboard events
    SoEventCallback* evCB = new SoEventCallback;
    evCB->addEventCallback(SoKeyboardEvent::getClassTypeId(), keyboardCB, nullptr);
    root->addChild(evCB);

    // Viewer
    SoQtExaminerViewer* viewer = new SoQtExaminerViewer(window);
    viewer->setSceneGraph(root);
    viewer->setTitle("3D Rope Simulation");
    viewer->setBackgroundColor(SbColor(0.93f, 0.93f, 0.97f));
    viewer->show();

    updateScene();
    cam->viewAll(root, viewer->getViewportRegion());

    // Focus the main window, not the viewer object
    window->setFocusPolicy(Qt::StrongFocus);
    window->setFocus();

    std::printf("\nControls:\n");
    std::printf("  Click once inside the 3D window first\n");
    std::printf("  W   -> toggle wind ON/OFF\n");
    std::printf("  R   -> reset rope\n");
    std::printf("  Esc -> quit\n\n");
    std::fflush(stdout);

    // Timer
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

    root->unref();
    return 0;
}