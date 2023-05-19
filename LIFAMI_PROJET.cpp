// Mon mini-projet combien maths et physique c'est incroyable non ? Il y aurat :
// - une bille de fliper et ses collisons
// - un ressort pour lancer la balle en début de partie
// - une horloge et ses aiguilles avec ses collisions
// - une magnifique image de dauphin
// J'ai trouvé de la documentation ici :
// - https://perso.liris.cnrs.fr/alexandre.meyer/grapic/html/index.html
// - http://fou-de-flippers.e-monsite.com/pages/les-flippers/petit-lexique-sur-les-organes-du-jeu.html pour le vocabulaire
// - https://lexique.netmath.ca/distance-entre-un-point-et-une-droite/ comme sont nom l'indique
// - https://pub.phyks.me/sdz/sdz/eorie-des-collisions.html#Thoriedescollisions : collisions

#include <Grapic.h>
#include <math.h>
#include <iostream>
using namespace grapic;
using namespace std;

const int DIMH = 600;
const int DIML = 400;
const float G = 35;
const float dt = 0.02;
const float FRICTION = 0.9;

//*************************************************************************************************//
//-------------------------OPERATEURS DE BASE & STRUCTURES DE DONNEES-----------------------------//
//***********************************************************************************************//

struct Vec2
{
    float x,y;
};

Vec2 make_Vec2 (float x, float y)
{
    Vec2 res;
    res.x = x;
    res.y = y;
    return res;
}

Vec2 make_Vec2_exp (float r ,float theta_deg)
{
    Vec2 res;
    res.x = r*cos (theta_deg*M_PI/180);
    res.y = r*sin (theta_deg*M_PI/180);
    return res;
}

Vec2 operator + (Vec2 opg , Vec2 opd)
{
    Vec2 res;
    res.x=opg.x + opd.x;
    res.y=opg.y + opd.y;
    return res;
}

Vec2 operator - (Vec2 opg , Vec2 opd)
{
    Vec2 res;
    res.x=opg.x - opd.x;
    res.y=opg.y - opd.y;
    return res;
}

Vec2 operator * (float lambda , Vec2 A)
{
    Vec2 res;
    res.x=lambda * A.x;
    res.y=lambda * A.y;
    return res;
}

Vec2 operator * (Vec2 A, float lambda)
{
    Vec2 res;
    res.x=A.x * lambda;
    res.y=A.y * lambda;
    return res;
}

Vec2 operator * (Vec2 opg, Vec2 opd)
{
    Vec2 res;
    res.x=opg.x*opd.x - opg.y*opd.y;
    res.y=opg.x*opd.y + opg.y*opd.x;
    return res;
}

struct bille
{
    Vec2 p,F,v; // la position, la force, la vitesse
    float m;  //la masse
    Image ImBille;
};

struct fond_obstacle
{
    Image ImHorloge;
    Image ImFond;

    Vec2 angle_b_g_1;
    Vec2 angle_b_g_2;

    Vec2 angle_b_d_1;
    Vec2 angle_b_d_2;

    Vec2 angle_h_g_1;
    Vec2 angle_h_g_2;

    Vec2 angle_h_d_1;
    Vec2 angle_h_d_2;

    Vec2 couloir1;
    Vec2 couloir2;

    Vec2 aiguille_h_1;
    Vec2 aiguille_h_2;

    Vec2 aiguille_m_1;
    Vec2 aiguille_m_2;

    Vec2 obstacle1;
    Vec2 obstacle2;
    Vec2 obstacle3;

};

struct ressort
{
    Vec2 pos;
    Image ImRessort;
};

struct batteur
{
    bool on_g;
    Vec2 p_g;
    Vec2 p_g2;
    float time_g;

    bool on_d;
    Vec2 p_d;
    Vec2 p_d2;
    float time_d;
};

struct aiguille
{
    Vec2 p_h1;
    Vec2 p_h2;
    float seq_h;

    Vec2 p_m1;
    Vec2 p_m2;
    float seq_m;
};

struct World
{
    fond_obstacle f;
    bille b;
    aiguille a;
    batteur bat;
    ressort r;

    int score;
    int score_fin;
    float time_start;
    float time_end;
    bool activ;
};

//*************************************************************************************************//
//--------------------------------PROCEDURE & FONCTION AVANCE-------------------------------------//
//***********************************************************************************************//

void partAddForce(bille &p, Vec2 f)
{
    p.F=p.F+f;
}

Vec2 rotation (Vec2 A, float cx, float cy, float R)
{
    Vec2 res;
    Vec2 t = make_Vec2(cx,cy);
    Vec2 rot = make_Vec2_exp(1,R);
    res = (A-t) * rot + t;
    return res;
}

float norme (Vec2 A)
{
    return sqrt(A.x*A.x+A.y*A.y);
}

float Distance_point (Vec2 A, Vec2 B) // distance entre deux points
{
    return sqrt( ((B.x-A.x)*(B.x-A.x)) + ((B.y-A.y)*(B.y-A.y)));
}

bool Distance_droit (Vec2 A, Vec2 B, Vec2 pos) // renvoie vraie si distance a pos est inférieur a 8
{
    Vec2 u; //ici on calul le vecteur de la doite AB
    u = make_Vec2(B.x - A.x, B.y - A.y);

    Vec2 AP; // permet de calculer la projection orthogonal de pos sur la droite
    AP = make_Vec2(pos.x - A.x, pos.y - A.y);

    float num = u.x*AP.y - u.y*AP.x; //norme du vecteur v (droite AB)

    if (num < 0)
    {
        num = -num; // si négatif on prend l'opposer
    }

    float deno = norme(u);
    float DISTANCE = num / deno; // donc la distance est la division de la norme de la projection orthogonal par norme de la doite AB

    if (DISTANCE<8) // si la distance est inférieur au rayon il y'a collision
    {
        return true;
    }
    else return false;
}

float prodscal (Vec2 u, Vec2 v)
{
    return u.x*v.x+u.y*v.y;
}

Vec2 Projection (Vec2 A, Vec2 B, Vec2 C) //projection de C sur la droite AB
{
    Vec2 AB, AC, res;
    AB = make_Vec2(B.x-A.x, B.y-A.y);
    AC = make_Vec2(C.x-A.x, C.y-A.y);

    float t = prodscal(AB,AC) / prodscal(AB,AB); // t est le coifficien de l'equation de la droite D(t) = A + t*AB;
    res = make_Vec2( (A.x + t*AB.x), (A.y + t*AB.y));
    return res;
}


bool CollisionPointCercle (Vec2 A, Vec2 pos) //renvoie vrais si la distance entre le point et le centre du cercle est inférieur au rayon
{
    int d = Distance_point(A,pos);

    if (d > 8)
    {
        return false;
    }
    else return true;
}

bool CollisionSegment (Vec2 A, Vec2 B, Vec2 pos) // renvoie vrais si la balle plus son rayon touche le segment AB
{
    if (Distance_droit(A,B,pos) == 0)
    {
        return false;
    }

    Vec2 AB,AP,BP;
    AB = make_Vec2(B.x-A.x, B.y-A.y);
    AP = make_Vec2(pos.x-A.x, pos.y-A.y);
    BP = make_Vec2(pos.x-B.x, pos.y-B.y);

    float prodscal1 = prodscal(AB,AP);
    float prodscal2 = prodscal(-1*AB,BP);

    if (prodscal1 >= 0 && prodscal2 >= 0) // on test si pos est entre A et B
    {
        return true;
    }

    if (CollisionPointCercle(A,pos)) // si A est dans le cercle de centre pos et de rayon 8
    {
        return true;
    }

    if (CollisionPointCercle(B,pos)) // si B est dans le cercle
    {
        return true;
    }
    else false;
}

Vec2 Normal (Vec2 A, Vec2 B, Vec2 pos) // revoie la normal normalisée (avec un norme = 1)
{
    Vec2 AB, AP, N;

    AB = make_Vec2(B.x-A.x, B.y-A.y);
    AP = make_Vec2(pos.x-A.x, pos.y-A.y);

    N.x = -AB.y*((AB.x*AP.y) - (AB.y*AP.x)); // en utilisant les deux produits vectoriels
    N.y = AB.x*((AB.x*AP.y) - (AB.y*AP.x));  // on trouve les formules suivantes

    N = (1.0 / norme(N))*N;

    return N;
}

Vec2 VecteurRebond (Vec2 u, Vec2 N) // renvoie le vecteur rebond obtenue grace au produit scalaire
{
    Vec2 res;
    float pdscal = prodscal(u,N); // nouveau vecteur proportionnelle au vecteur d'arriver u
    res = make_Vec2((u.x -2*pdscal*N.x),(u.y -2*pdscal*N.y));
    return res;
}

//*************************************************************************************************//
//---------------------------------INITIALISATION & AFFICHAGE-------------------------------------//
//***********************************************************************************************//

void init_bille (bille &b)
{
    b.m = 0.080; // masse 80 grammes
    b.p = make_Vec2(381,80); //(378 mid ressort,80)
    b.v = make_Vec2(0,0);
    b.ImBille = image("data/memory/bille1.png");
}

void init_fond (fond_obstacle &fo)
{
    fo.ImHorloge = image("data/memory/horloge2.png");
    fo.ImFond = image("data/memory/fond.jpg");

    fo.angle_b_d_1 = make_Vec2(360,125); // initialisation des deux point
    fo.angle_b_d_2 = make_Vec2(260,25);

    fo.angle_b_g_1 = make_Vec2(25,125); //Bas gauche
    fo.angle_b_g_2 = make_Vec2(125,25);

    fo.angle_h_d_1 = make_Vec2(DIML,500); //Haut droit
    fo.angle_h_d_2 = make_Vec2(300,DIMH);

    fo.angle_h_g_1 = make_Vec2(0,500); //Haut gauche
    fo.angle_h_g_2 = make_Vec2(100,DIMH);

    fo.couloir1 = make_Vec2(359,0);
    fo.couloir2 = make_Vec2(359,500);

    fo.obstacle1 = make_Vec2(75,450);
    fo.obstacle2 = make_Vec2(270,425);
    fo.obstacle3 = make_Vec2(200,550);

}

void init_aiguille (aiguille &a)
{
    a.p_h1=make_Vec2(DIML/2,DIMH/2);
    a.p_h2=make_Vec2((DIML/2)+50,(DIMH/2)+50);
    a.seq_h = 4;

    a.p_m1=make_Vec2(DIML/2,DIMH/2);
    a.p_m2=make_Vec2((DIML/2)-50, (DIMH/2)+75);
    a.seq_m = 1;
}

void init_ressort (ressort &r)
{
    r.ImRessort = image("data/memory/ressort.jpg");
    r.pos = make_Vec2(360,0);
}

void init_batteur (batteur &b)
{
    b.on_g = 0;
    b.p_g = make_Vec2(125,25);
    b.p_g2 = make_Vec2(170,5);
    b.time_g =0;

    b.on_d = 0;
    b.p_d = make_Vec2(260,25);
    b.p_d2 = make_Vec2(215,5);
    b.time_d = 0;
}

void init_World (World &w)
{
    w.score = 0;
    w.time_start = 0;
    w.time_end = 0;
    w.activ = 1;

    init_bille(w.b);
    init_fond(w.f);
    init_ressort(w.r);
    init_aiguille(w.a);
    init_batteur(w.bat);
}

void draw_bille (bille b)
{
    image_draw(b.ImBille, b.p.x-8, b.p.y-8, 16 , 16);
}

void draw_batteur (batteur b)
{
    if (b.on_d == 0)
    {
        color(42,62,164);
        line(b.p_d.x, b.p_d.y, b.p_d2.x, b.p_d2.y);
    }
    else
    {
        b.p_d2 = rotation(b.p_d2, b.p_d.x,b.p_d.y, -40);

        color(42,62,164);
        line(b.p_d.x, b.p_d.y, b.p_d2.x, b.p_d2.y);
    }

    if (b.on_g == 0)
    {
        color(42,62,164);
        line(b.p_g.x,b.p_g.y, b.p_g2.x,b.p_g2.y);
    }
    else
    {
        b.p_g2 = rotation(b.p_g2, b.p_g.x,b.p_g.y, +40);
        color(42,62,164);
        line(b.p_g.x, b.p_g.y, b.p_g2.x, b.p_g2.y);
    }
}

void draw_fond(fond_obstacle fo)
{
    image_draw(fo.ImFond, (DIML/2)-55, (DIMH/2)-75 , 120, 145);
    image_draw(fo.ImHorloge, (DIML/2)-115, (DIMH/2)-115, 230, 230);


    color(42,62,164); // couloir gauche;
    line(fo.couloir1.x,fo.couloir1.y,fo.couloir2.x,fo.couloir2.y);

    color(34,49,126); //triangle bas gauche
    triangleFill(25,125, 125,25, 25,25);

    color (34,49,126); //triangle bas droite
    triangleFill(360,125, 260,25, 360,25);

    color(42,62,164); //triangle haut gauche
    triangleFill(0,500, 100,DIMH, 0,DIMH);

    color(42,62,164); //triangle haut droit
    triangleFill(DIML,500, 300,DIMH, DIML,DIMH);

    color(34,49,126); //obstacle1
    triangleFill(fo.obstacle1.x-15, fo.obstacle1.y-15, fo.obstacle1.x+15, fo.obstacle1.y-15, fo.obstacle1.x, fo.obstacle1.y+15);

    color(34,49,126); //obstacle2
    triangleFill(fo.obstacle2.x-15, fo.obstacle2.y-15, fo.obstacle2.x+15, fo.obstacle2.y-15, fo.obstacle2.x, fo.obstacle2.y+15);

    color(34,49,126); //obstacle3
    triangleFill(fo.obstacle3.x-15, fo.obstacle3.y-15, fo.obstacle3.x+15, fo.obstacle3.y-15, fo.obstacle3.x, fo.obstacle3.y+15);
}

void draw_deco ()
{
    color(255,255,255);
    triangleFill(0,550, 50,DIMH, 0,DIMH);

    color(255,255,255);
    triangleFill(50,100, 100,50, 50,50);
    color(34,49,126);
    line(125,25, 25,125);
    color(34,49,126);
    line(120,25, 25,120);

    color(255,255,255);
    triangleFill(335,100, 275,50, 335,50);
    color(34,49,126);
    line(260,25, 360,125);
    color(34,49,126);
    line(265,25, 360,120);

    color(255,255,255);
    circleFill(75,445,7);

    color(255,255,255);
    circleFill(270,420,7);

    color(255,255,255);
    circleFill(200,545,7);

    color(42,62,164); //ligne gauche
    line(0,0,0,DIMH);

    color(42,62,164); //ligne droit
    line(DIML-1,0,DIML-1,DIMH);

     color(42,62,164); //ligne haut
    line(0,DIMH-1,DIML,DIMH-1);

}

void draw_aiguille (aiguille a)
{
    color(0,0,0);
    line(a.p_h1.x, a.p_h1.y, a.p_h2.x, a.p_h2.y);

    color(0,0,0);
    line(a.p_m1.x, a.p_m1.y, a.p_m2.x, a.p_m2.y);
}

void draw_ressort (ressort r)
{
    image_draw(r.ImRessort, r.pos.x, r.pos.y, 40, 50);
}

void draw_game_over (World w)
{
    color(0,0,0);
    rectangleFill(0,0,DIML,DIMH);

    color(236,15,290);
    print((DIML/2)-50,DIMH/2,"GAME OVER");

    color(255,255,255);
    print((DIML/2)-90,(DIMH/2)-20,"Temps de jeu :");
    color(255,255,255);
    print((DIML/2)+50,(DIMH/2)-20,w.time_end);

    color(255,255,255);
    print((DIML/2)-140,(DIMH/2)-120,"Appuyez sur ESPACE pour rejouer");

    color(255,255,255);
    print(288,550,"Score : ");
    color(255,255,255);
    print(365, 550, w.score_fin);
}

void draw_World (World w, bool activer)
{
    if (activer == 1)
    {
        draw_ressort(w.r);
        draw_batteur(w.bat);
        draw_fond(w.f);
        draw_aiguille(w.a);
        draw_deco();
        draw_bille(w.b);
    }
    else draw_game_over(w);
}


//*************************************************************************************************//
//----------------------------------UPDATE & COLLISIONS-------------------------------------------//
//***********************************************************************************************//


void update_bille (bille &b)
{
    b.F.x = 0;
    b.F.y = 0;

    Vec2 force = make_Vec2(0, -G*b.m); // ici on crée un vecteur qui joue le role de la gravité sur terre.
    partAddForce(b, force);
    b.v = b.v+1/b.m*dt*b.F; // formule de la vitesse en fonction du temps dt, de la masse b.m et de la force b.F.
    b.p = b.p + dt * b.v; // formule de la position en fonction du temps et du vecteur vitesse.


    if (b.p.y < 50 && b.p.x > 360) //collision pour le ressort
    {
        b.p.y=50;
        b.v.y=-b.v.y;
        b.v = b.v*FRICTION;
    }

    if (b.p.x <8)
    {
        b.p.x=b.p.x;
        b.v.x=-b.v.x;
    }

    if (b.p.x > DIML-8)
    {
        b.p.x=b.p.x;
        b.v.x=-b.v.x;
    }

    if (b.p.y > DIMH-8)
    {
        b.p.y=b.p.y;
        b.v.y=-b.v.y;
    }

    if (isKeyPressed(SDLK_UP) && b.p.x > 360 && b.p.y < 60) // lancement de la bille
    {
        b.v=make_Vec2(0,210);
    }

    Vec2 droite1 = make_Vec2(359,0);
    Vec2 droite2 = make_Vec2(359,450);
    if (CollisionSegment(droite1,droite2,b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(droite1, droite2, b.p));
        b.v.y = b.v.y * FRICTION;
    }
}

void update_fond (fond_obstacle f, bille &b, World &w)
{
    if (CollisionSegment(f.angle_h_g_1,f.angle_h_g_2, b.p) == 1) //collision angle haut gauche
    {
        b.v = VecteurRebond(b.v, Normal(f.angle_h_g_1, f.angle_h_g_2, b.p));
    }

    if (CollisionSegment(f.angle_h_d_1, f.angle_h_d_2, b.p) == 1) //collison angle haut droit
    {
        b.v = VecteurRebond(b.v,Normal(f.angle_h_d_1, f.angle_h_d_2, b.p));
    }

    if (CollisionSegment(f.angle_b_d_1, f.angle_b_d_2, b.p) == 1) //collision angle bas droit
    {
        b.v = VecteurRebond(b.v,Normal(f.angle_b_d_1, f.angle_b_d_2, b.p));
    }

    if (CollisionSegment(f.angle_b_g_1, f.angle_b_g_2, b.p) == 1) //collision angle bas gauche
    {
        b.v = VecteurRebond(b.v,Normal(f.angle_b_g_1, f.angle_b_g_2, b.p));
    }

    Vec2 bas_g = make_Vec2(f.angle_b_g_1.x, f.angle_b_g_1.y - 100); // Collision coté du triangle bas gauche
    if (CollisionSegment(f.angle_b_g_1, bas_g, b.p))
    {
        b.v = VecteurRebond(b.v, Normal(f.angle_b_g_1, bas_g, b.p));
    }

    //Collision obstacle 1
    Vec2 obs1_bg = make_Vec2(f.obstacle1.x -15, f.obstacle1.y -15);
    Vec2 obs1_bd = make_Vec2(f.obstacle1.x +15, f.obstacle1.y -15);
    Vec2 obs1_h = make_Vec2(f.obstacle1.x, f.obstacle1.y +15);

    if (CollisionSegment(obs1_bd,obs1_bg, b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(obs1_bd, obs1_bg, b.p));
        w.score = w.score + 10;
    }
    if (CollisionSegment(obs1_bd,obs1_h, b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(obs1_bd, obs1_h, b.p));
        w.score = w.score + 10;
    }
    if (CollisionSegment(obs1_bg,obs1_h, b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(obs1_bg, obs1_h, b.p));
        w.score = w.score + 10;
    }

    //Collsion obstacele 2
    Vec2 obs2_bg = make_Vec2(f.obstacle2.x -15, f.obstacle2.y -15);
    Vec2 obs2_bd = make_Vec2(f.obstacle2.x +15, f.obstacle2.y -15);
    Vec2 obs2_h = make_Vec2(f.obstacle2.x, f.obstacle2.y +15);

    if (CollisionSegment(obs2_bd,obs2_bg, b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(obs2_bd, obs2_bg, b.p));
        w.score = w.score + 10;
    }
    if (CollisionSegment(obs2_bd,obs2_h, b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(obs2_bd, obs2_h, b.p));
        w.score = w.score + 10;
    }
    if (CollisionSegment(obs2_bg,obs2_h, b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(obs2_bg, obs2_h, b.p));
        w.score = w.score + 10;
    }

    //Collsion obstacele 3
    Vec2 obs3_bg = make_Vec2(f.obstacle3.x -15, f.obstacle3.y -15);
    Vec2 obs3_bd = make_Vec2(f.obstacle3.x +15, f.obstacle3.y -15);
    Vec2 obs3_h = make_Vec2(f.obstacle3.x, f.obstacle3.y +15);

    if (CollisionSegment(obs3_bd,obs3_bg, b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(obs3_bd, obs3_bg, b.p));
        w.score = w.score + 10;
    }
    if (CollisionSegment(obs3_bd,obs3_h, b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(obs3_bd, obs3_h, b.p));
        w.score = w.score + 10;
    }
    if (CollisionSegment(obs3_bg,obs3_h, b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(obs3_bg, obs3_h, b.p));
        w.score = w.score + 10;
    }
}

void update_batteur (batteur &bat, bille &b, World w)
{
    if (isKeyPressed(SDLK_LEFT)) //gauche
    {
        bat.on_g = 1;
        bat.time_g = w.time_start;
    }

    if (bat.on_g == 1 && w.time_start-bat.time_g > 0.3) // si le batteur gauche est au repos on applique telle collision
    {
        bat.on_g = 0;
        bat.time_g = 0;
    }
    // Collision gauche
    if (bat.on_g == 0 && CollisionSegment(bat.p_g, bat.p_g2, b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(bat.p_g, bat.p_g2, b.v));
    }

    if (bat.on_g == 1 && CollisionSegment(bat.p_g, bat.p_g2, b.p) == 1)
    {
        Vec2 point_g2 = rotation(bat.p_g2, bat.p_g.x,bat.p_g.y, +40);
        b.v = VecteurRebond(b.v, Normal(bat.p_g, point_g2, b.v));

    }

    if (isKeyPressed(SDLK_RIGHT)) //droite
    {
        bat.on_d = 1;
        bat.time_d = w.time_start;
    }

    if (bat.on_d == 1 && w.time_start-bat.time_d > 0.3)
    {
        bat.on_d = 0;
        bat.time_d = 0;
    }

    // Collision droite
    if (bat.on_d == 0 && CollisionSegment(bat.p_d, bat.p_d2, b.p) == 1)
    {
        b.v = VecteurRebond(b.v, Normal(bat.p_d, bat.p_d2, b.v));

    }

    if (bat.on_d == 1 && CollisionSegment(bat.p_d, bat.p_d2, b.p) == 1)
    {
        Vec2 point_d2 = rotation(bat.p_d2, bat.p_d.x,bat.p_d.y, +40);
        b.v = VecteurRebond(b.v, Normal(bat.p_g, point_d2, b.v));
    }
}

void update_aiguille (aiguille &a, bille &b, float time)
{

    if  (time >= a.seq_h)
    {
        a.p_h2 = rotation(a.p_h2, a.p_h1.x, a.p_h1.y, -45);
        a.seq_h = a.seq_h + 4;
    }
    if (CollisionSegment(a.p_h1, a.p_h2, b.p) == 1)
    {
        b.v = VecteurRebond(b.v,Normal(a.p_h1, a.p_h2, b.p));
    }

    if  (time >= a.seq_m)
    {
        a.p_m2 = rotation(a.p_m2, a.p_m1.x, a.p_m1.y, -20);
        a.seq_m = a.seq_m + 1;
    }

    if (CollisionSegment(a.p_m1, a.p_m2, b.p) == 1)
    {
        b.v = VecteurRebond(b.v,Normal(a.p_m1, a.p_m2, b.p));
    }
}

void update_World (World &w, bool activer)
{
    if (activer == 1)
    {
        color(255,255,255);
        print(365, 550, w.score);

        w.time_start = elapsedTime();

        update_bille(w.b);
        update_fond(w.f, w.b, w);
        update_batteur(w.bat, w.b, w);
        update_aiguille(w.a, w.b, w.time_start);
    }
    else
    {
        w.score_fin = w.score;
        w.time_end = w.time_start;
    }

}

int main(int , char** )
{
    bool stop=false;

    winInit("Flip'heure", DIML, DIMH);

    backgroundColor( 255, 255, 255 );
    World w;
    float t = 15;
    init_World(w);

    while( !stop )
    {
        winClear();
        draw_World(w,w.activ);
        update_World(w,w.activ);

        if (w.time_start > t && w.b.p.x < 360 && w.b.p.y > 50) //ajoute 100 points toute les 15 seconde
        {
            w.score = w.score + 100;
            t = t+20;
        }

        if (w.b.p.y < 0)
        {
            w.activ = 0;
        }

        if (isKeyPressed(SDLK_SPACE) && w.activ == 0) //replay
        {
            init_World(w);
            w.time_start = 0;
            float t = 15;
        }

        stop = winDisplay();
    }
    winQuit();
    return 0;
}
