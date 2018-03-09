// Determine the projection of a point on a plane centered at c1 with normal n1
function project(p1, c1, n1) =
    p1 - (n1 * (p1 - c1)) * n1 / (n1 * n1);

// Determine the angle between two points and a centre in 3D space    
// c^2 = a^2 + b^2 -2ab * cos(C)
// <=> cos(C) = (a^2 + b^2 -c^2) / (2ab)
function getAngle(p1, c1, p2) =
    acos(((p1-c1)*(p1-c1) + (p2-c1)*(p2-c1) - (p1-p2)*(p1-p2)) / 
        (2*norm(p1-c1)*norm(p2-c1)));

// Generate a line between two points in 3D space
module line3D(tp1,tp2, tk=1, dp=1){
    p1 = c3D(tp1);
    p2 = c3D(tp2);
    sRot = (rToS(p1-p2));
    ll = norm(p1-p2);
    translate((p1+p2)/2)
        rotate(sRot) rotate([90,0,0]) rotate([90,90,0])
            cylinder(h=ll, d=dp, center=true);
}

// convert a rotation angle to a rotation matrix
function rot2Mat(rotVec, axis) =
    (len(rotVec) == 2) ?
        rot2Mat([rotVec[0], rotVec[1], 0], axis) :
    (axis == "x") ?
        [[1,              0,               0],
         [0, cos(rotVec[0]),  sin(rotVec[0])],
         [0, sin(rotVec[0]), -cos(rotVec[0])]] :
    (axis == "y") ?
        [[ cos(rotVec[1]), 0, sin(rotVec[1])],
         [              0, 1,              0],
         [-sin(rotVec[1]), 0, cos(rotVec[1])]] :
    (axis == "z") ?
        [[ cos(rotVec[2]), sin(rotVec[2]), 0],
         [-sin(rotVec[2]), cos(rotVec[2]), 0],
         [0,              0,               1]] : undef;

// rotate a point (or points)
function myRotate(rotation, points) =
    (len(points[0]) == undef) ?
        myRotate(rotation, [points])[0] :
        c3D(points) * rot2Mat(rotation, "x")
               * rot2Mat(rotation, "y")
               * rot2Mat(rotation, "z");

// translate a point (or points)
function myTranslate(ofs, points, acc = []) =
    (len(points[0]) == undef) ?
        myTranslate(ofs, [points])[0] :
        (len(acc) == len(points)) ? acc :
            myTranslate(ofs, points, concat(acc, [points[len(acc)] + ofs]));

// convert point to 3D by setting Z to zero (if not present)
function c3D(tPoints) = 
    (len(tPoints[0]) == undef) ?
        c3D([tPoints])[0] :
        (len(tPoints[0]) < 3) ?
            tPoints * [[1,0,0],[0,1,0]] :
            tPoints;

// Determine spherical rotation for cartesian coordinates
function rToS(pt) = 
    [-acos((pt[2]) / norm(pt)), 
         0,
         -atan2(pt[0],pt[1])];

// Rotate a position around an angle
function rotPos(r, ang) =
    [r * cos(ang), r * sin(ang), 0];


function intp(p1, p2, thr=0.5, res = []) = 
    (norm(p2-p1) <= thr) ? concat(res,[p1]) :
        intp(p1=p1 + (thr/norm(p2-p1)) * (p2-p1), p2=p2, 
            thr=thr, res = concat(res,[p1]));

// see https://stackoverflow.com/questions/14066933/
//        direct-way-of-computing-clockwise-angle-between-2-vectors
//    dot = p1 * p2;
//    det = (p1[0]*p2[1]*n1[2] + p2[0]*n1[1]*p1[2] + n1[0]*p1[1]*p2[2]) -
//          (n1[0]*p2[1]*p1[2] + p1[0]*n1[1]*p2[2] + p2[0]*p1[1]*n1[2]);
//    atan2(det, dot);

// determine angle between two points with a given normal orientation
function getNormAngle(p1, n1, p2) =
    atan2((p1[0]*p2[1]*n1[2] + p2[0]*n1[1]*p1[2] + n1[0]*p1[1]*p2[2]) -
          (n1[0]*p2[1]*p1[2] + p1[0]*n1[1]*p2[2] + p2[0]*p1[1]*n1[2]), p1 * p2);

// determine angle between two points and a centre with a given normal orientation
function getNPAngle(p1, c1, n1, p2) = 
    getNormAngle(p1=p1-c1, n1=n1 / norm(n1), p2=p2-c1);

// calculate offset based on a given array length, wrapping around to zeroth element
function arrMod(arrBig, arrSmall, ofs) =
    arrBig[(len(arrSmall) + len(arrBig) + ofs) % len(arrBig)];

// calculate offset based on a given length, wrapping around to zeroth element
function wrapMod(bigLength, arrLength, ofs) =
    (arrLength + bigLength + ofs) % bigLength;

//        t0 = p0 + myRotate(rToS(rPlanes[0]), myRotate([0,0,-rawPreRots[0]], c3D(myPoints[0])));
//        tm1 = pm1 + myRotate(rToS(rPlanes[len(rPlanes)-1]), 
//                myRotate([0,0,-rawPreRots[len(rawPreRots)-1]], c3D(myPoints[0])));
//        pt0 = project(p1=t0, c1=pm1, n1=(pm1-p0));
//        lfAng = -getNPAngle(p1 = pt0, c1 = pm1, n1=(pm1-pm2), p2=tm1);


// work out planar rotations for path slices, minimising distance between the
// first coordinate in the polygon
function getPreRotations(extrudePath, refPt, polyNormals, merge=false, prs = [0]) =
    (len(prs) >= (len(extrudePath))) ? prs :
        getPreRotations(extrudePath=extrudePath, refPt=refPt, 
            polyNormals=polyNormals, merge=merge,
            prs=concat(prs,getNPAngle(p1=project(
                p1=arrMod(extrudePath,prs,-1) + 
                    myRotate(rToS(arrMod(polyNormals,prs,-1)),
                    myRotate([0,0,-prs[len(prs)-1]], refPt)),
                n1=arrMod(polyNormals,prs,-1),
                c1=arrMod(extrudePath, prs, -1)),
            c1=arrMod(extrudePath,prs, 0),
            n1=arrMod(polyNormals,prs, 0),
            p2=arrMod(extrudePath,prs, 0) + 
                myRotate(rToS(arrMod(polyNormals,prs,0)), refPt))));
                    
// spreads an adjustment across all values in an array to reduce jumps
function spreadError(a, adj, acc = []) = 
   (len(acc) == len(a)) ? acc :
        spreadError(a = a, adj = adj,
            acc = concat(acc, a[len(acc)] + (adj / (len(a))) * len(acc)));

function getRotationNormals(polyPath, merge = false, acc=[], aDone = 0) =
    (aDone >= len(polyPath)) ? acc :
        getRotationNormals(polyPath = polyPath, merge=merge, 
            acc=concat(acc, [arrMod(polyPath, acc, ((!merge) && (aDone>=(len(polyPath)-1)))?0:1) - 
                arrMod(polyPath,acc,((!merge) && (aDone==0))?0:-1)]),
            aDone=aDone+1);

// set up massive point array for polyhedron
function makePolyPoints(polyPath, polyForm, polyAngles, polyNormals, merge = false, acc = [], aDone = 0) =
  (aDone >= len(polyPath)) ? acc :
    makePolyPoints(polyPath=polyPath, polyForm=polyForm, polyAngles=polyAngles, polyNormals=polyNormals,
        merge = merge,
        acc=concat(acc,  
            [myTranslate(arrMod(polyPath,acc,0),
                myRotate(rToS(polyNormals[aDone]),
                myRotate([0,0,-polyAngles[len(acc)]], c3D(polyForm))))]),
        aDone=aDone + 1);

// removes the top level array from the array A
function flatten(A, acc = [], aDone = 0) =
    (aDone >= len(A)) ? acc :
       flatten(A, acc=concat(acc, A[aDone]), aDone = aDone + 1);

// creates a triangle list joining adjacent polygons
function makeTriAdjs(pathLen, formLen, i, acc = [], aDone = 0) =
    (aDone >= formLen) ? acc :
        makeTriAdjs(pathLen, formLen, i, 
            acc = concat(acc, [[
              [(i*formLen + wrapMod(formLen, aDone, 1)) % (pathLen*formLen),
              (i*formLen + wrapMod(formLen, aDone, 0)) % (pathLen*formLen), 
              (i*formLen + wrapMod(formLen, aDone, 1) + formLen) % (pathLen*formLen)],
              [(i*formLen + wrapMod(formLen, aDone, 1) + formLen) % (pathLen*formLen),
              (i*formLen + wrapMod(formLen, aDone, 0)) % (pathLen*formLen), 
              (i*formLen + wrapMod(formLen, aDone, 0) + formLen) % (pathLen*formLen)]]]),
            aDone = aDone+1);

myPath = [ for(t = [0:(360 / 51):359]) [ // trefoil knot
    5*(.41*cos(t) - .18*sin(t) - .83*cos(2*t) - .83*sin(2*t) - .11*cos(3*t) + .27*sin(3*t)),
    5*(.36*cos(t) + .27*sin(t) - 1.13*cos(2*t) + .30*sin(2*t) + .11*cos(3*t) - .27*sin(3*t)),
    5*(.45*sin(t) - .30*cos(2*t) +1.13*sin(2*t) - .11*cos(3*t) + .27*sin(3*t))] ];

//myPath = [[-1,0,0],[1,0,0],[2,1,0.5],[2,3,1.5],
//    [1,4,2],[-1,4,3],[-2,3,3.5],[-2,1,4.5],[-1,0,5]]; // pentagon spiral

//myPath = [ for(t = [0:(360/20):359]) 8 * [cos(t+45),sin(t+45),0] ]; // octagon

ofs1=15;
    
myPoints = [ for(t = [0:(360/8):359]) ((t==90)?1:2) * [cos(t+ofs1),sin(t+ofs1)]];
myPointsChunk = [ for(t = [45:(360/8):136]) ((t==90)?1.5:1.9) * [cos(t+ofs1),sin(t+ofs1)]];
//myPoints = [ for(t = [0:(360/8):359]) 2 * [cos(t+45),sin(t+45)]];

module path_extrude(exPath = myPath, exShape = myPoints, exRots = [0], merge=false){
    if((exShape == undef) || (exPath == undef)){
        echo("Path or shape not defined");
    } else {
        rPlanes = getRotationNormals(exPath, merge=merge);
        // calculate rotations to reorient polygons to best match consecutive copies
        rawPreRots = getPreRotations(extrudePath=exPath, refPt=c3D(exShape[0]),
                   polyNormals=rPlanes, merge=merge, prs=exRots);
        // calculate rotation between last polygon and first polygon
        pp1 = exPath[1];
        p0 = exPath[0];
        pm1 = exPath[len(exPath)-1];
        pm2 = exPath[len(exPath)-2];
        t0 = p0 + myRotate(rToS(rPlanes[0]), myRotate([0,0,-rawPreRots[0]], c3D(exShape[0])));
        tm1 = pm1 + myRotate(rToS(rPlanes[len(rPlanes)-1]), 
                myRotate([0,0,-rawPreRots[len(rawPreRots)-1]], c3D(exShape[0])));
        pt0 = project(p1=t0, c1=pm1, n1=rPlanes[len(rPlanes)-1]);
        lfAng = -getNPAngle(p1 = pt0, c1 = pm1, n1=rPlanes[len(rPlanes)-2], p2=tm1);
        preRots = (merge) ? spreadError(rawPreRots, -lfAng) : rawPreRots;
        polyPoints = flatten(makePolyPoints(polyPath=exPath, polyForm=exShape, 
            polyAngles=preRots, polyNormals=rPlanes, merge=merge));
        if(merge){
            polyhedron(points = polyPoints, 
                faces = flatten([ for(i = [0:(len(exPath)-1)]) 
                    flatten(makeTriAdjs(len(exPath), len(exShape), i)) ]));
        } else {
            polyhedron(points = polyPoints, 
                faces = concat(
                    concat(flatten([ for(i = [0:(len(exPath)-2)]) 
                        flatten(makeTriAdjs(len(exPath), len(exShape), i)) ]),
                        [[for(i = [0:(len(exShape)-1)]) i]]),
                    [[for(i = [0:(len(exShape)-1)]) ((len(exPath)*len(exShape))-1-i)]]));
        }
    }
}

path_extrude(exRots = [$t*360], merge=false);

color("lightblue") path_extrude(exRots = [$t*360], exShape=myPointsChunk, merge=false);