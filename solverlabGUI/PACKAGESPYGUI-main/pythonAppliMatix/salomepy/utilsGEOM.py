#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import math
import salomepy.xsalomesession as XSS

def getDXDYDZVector(vecteur):
    geompy = XSS.getGeompy()
    p1 = geompy.MakeVertexOnCurve(vecteur,0.)
    p2 = geompy.MakeVertexOnCurve(vecteur,1.)
    c1 = geompy.PointCoordinates(p1)
    c2 = geompy.PointCoordinates(p2)
    
    dx = c2[0]-c1[0]
    dy = c2[1]-c1[1]
    dz = c2[2]-c1[2]
    return dx, dy, dz
    
def OpposedVectors(uGeom, vGeom):
    eps = 0.00001
    u0, u1, u2 = getDXDYDZVector(uGeom)
    v0, v1, v2 = getDXDYDZVector(vGeom)
    ve0, ve1, ve2 = [u1*v2-u2*v1, u2*v0-u0*v2, u0*v1-u1*v0] #vectoriel
    scalaire = u0*v0 + u1*v1 + u2*v2
    normeVectoriel = math.sqrt(ve0*ve0 + ve1*ve1 + ve2*ve2)
    if normeVectoriel < eps and scalaire < 0: return True
    return False

def NormeVectorielle(uGeom):
    u0, u1, u2 = self.getDXDYDZVector(uGeom)
    return math.sqrt(u0*u0 + u1*u1 + u2*u2)

def ProduitScalaire(uGeom, vGeom):
    u0, u1, u2 = getDXDYDZVector(uGeom)
    v0, v1, v2 = getDXDYDZVector(vGeom)
    return u0*v0 + u1*v1 + u2*v2
    
def ProduitVectoriel(uGeom, vGeom):
    geompy = XSS.getGeompy()
    u0, u1, u2 = getDXDYDZVector(uGeom)
    v0, v1, v2 = getDXDYDZVector(vGeom)
    return geompy.MakeVectorDXDYDZ(u1*v2 - u2*v1, u2*v0 - u0*v2, u0*v1 - u1*v0)

def DiffLimitBetweenTwoVectors(u, v, seuilAngle, eps):
    if math.fabs(u[0]-v[0]) > seuilAngle:
        return False
    if math.fabs(u[1]-v[1]) > seuilAngle:
        return False
    if math.fabs(u[2]-v[2]) > seuilAngle:
        return False
    return True




######################################################
#from NUMODIS GenericFunctions.py
######################################################

def planeEquation( listvertex ):
    geompy = XSS.getGeompy()
    pt0=geompy.PointCoordinates(listvertex[0])
    pt1=geompy.PointCoordinates(listvertex[1])
    pt2=geompy.PointCoordinates(listvertex[2])
    vect01=[(pt1[0]-pt0[0]),(pt1[1]-pt0[1]),(pt1[2]-pt0[2])]
    vect02=[(pt2[0]-pt0[0]),(pt2[1]-pt0[1]),(pt2[2]-pt0[2])]
    vect=CrossProduct(vect01,vect02)
    nrm2= dNorm(vect)
    a=(vect[0]/nrm2) 
    b=(vect[1]/nrm2) 
    c=(vect[2]/nrm2)
    d=((pt0[0]*a)+(pt0[1]*b)+(pt0[2]*c))
    return [a, b, c, d]

def isInteger(s):
    try :
       i = int(s)
       return True
    except :
       return False

def calculateMinDistanceToFaces( vertex, facesList ):
    geompy = XSS.getGeompy()
    minDistances = []
    for i in range( len(facesList) ) :
       minDistances.append( geompy.MinDistance( vertex, facesList[i] ) )
    minDistances.sort()
    return minDistances[0]

def areParallel( vx,vy,vz, ux,uy,uz ):
    """This method returns True if the vectors v(vx,vy,vz) and u(ux,uy,uz) are parallel, False otherwise"""
    epsilon = 10e-6
    r1 = vx/ux
    r2 = vy/uy
    r3 = vz/uz
    if math.fabs(r1-r2) < epsilon and \
       math.fabs(r1-r3) < epsilon and \
       math.fabs(r2-r3) < epsilon:
      return True
    return False

def arePerpendicular( vx,vy,vz, ux,uy,uz ):
    """
    This method returns True if the vectors 
    v(vx,vy,vz) and u(ux,uy,uz) are perpendicular, 
    False otherwise
    """
    epsilon = 10e-6
    return vx*ux + vy*uy + vz*uz < epsilon

def haveTheSameDirection( u, v ):
    """
    This method returns True if vectors u and v 
    have the same direction, False otherwise.
    The given vectors must be colinear
    """
    geompy = XSS.getGeompy()
    Partition = geompy.MakePartition( [u], [v], [], [], geompy.ShapeType["VERTEX"] )
    interPoints = geompy.SubShapeAll( Partition, geompy.ShapeType["VERTEX"] )
    return len( interPoints ) == 2

def isInside( geomObjectA, geomObjectB ):
    """
    This method returns True if geomObjectA is inside geomObjectB,
    False otherwise
    """
    geompy = XSS.getGeompy()
    epsilon = 10e-6
    vertices = geompy.SubShapeAll( geomObjectA, geompy.ShapeType['VERTEX'] )
    for i in range( len(vertices) ):
       vertex = vertices[i]
       minDistance = geompy.MinDistance( geomObjectB, vertex ) 
       if minDistance > epsilon : return False
    return True

def getEquation( faceGeomObject ):
    """
    This method returns the equation of faceGeomObject as a dico
    of four pamars a, b, c and d
    """
    geompy = XSS.getGeompy()
    points = geompy.SubShapeAll( faceGeomObject, geompy.ShapeType["VERTEX"] )
    point_1, point_2, point_3 = points[0:3]

    point_1Coord = geompy.PointCoordinates( point_1 )
    point_2Coord = geompy.PointCoordinates( point_2 )
    point_3Coord = geompy.PointCoordinates( point_3 )

    u = geompy.MakeVector( point_1, point_2 )
    v = geompy.MakeVector( point_1, point_3 )
    u1, u2, u3 = getDXDYDZVector( u )
    v1, v2, v3 = getDXDYDZVector( v )

    normal = geompy.MakeVectorDXDYDZ( u2*v3-u3*v2, u3*v1-u1*v3, u1*v2-u2*v1 )
    #Calculating the equation of the plane defined by normal and containing point_1, point_2 and point_3
    a, b, c = getDXDYDZVector( normal )

    #normal must have the same direction as the exiting normal to the face given by geompy
    normalToFace = geompy.GetNormal( faceGeomObject, point_1 )
    if not OpposedVectors( normal, normalToFace ) :
       d = a*point_1Coord[0] + b*point_1Coord[1] + c*point_1Coord[2] 
       equation = { 'a':a, 'b':b, 'c':c, 'd':d }
    else : 
       d = -a*point_1Coord[0] - b*point_1Coord[1] - c*point_1Coord[2] 
       equation = { 'a':-a, 'b':-b, 'c':-c, 'd':d }
    return equation

def getContinuousListOfVertices( faceGeomObject ):
    """
    This method returns the list of the vertices of faceGeomObject 
    in a continuous way
    (i.e the face is parsed and the vertices (vertex) are added)
    """
    geompy = XSS.getGeompy()
    vertices = geompy.SubShapeAllIDs( faceGeomObject, geompy.ShapeType['VERTEX'] )
    edges = geompy.SubShapeAllIDs( faceGeomObject, geompy.ShapeType['EDGE'] )

    #Dico of correspondances edges-->vertices
    ES = {}
    for i in range( len(edges) ):
       edge = geompy.GetSubShape( faceGeomObject, [edges[i]] )
       extremities = []
       for j in range( len(vertices) ):
          vertex = geompy.GetSubShape( faceGeomObject, [vertices[j]] )
          if geompy.MinDistance( edge, vertex ) == 0 : extremities.append( geompy.GetSubShapeID(faceGeomObject, vertex) )
       ES[edges[i]] = extremities


    #Dico of correspondances vertices-->edges
    SE = {}
    for i in range( len(vertices) ):
       vertex = geompy.GetSubShape( faceGeomObject, [vertices[i]] )
       relatedEdges = []
       for j in range( len(edges) ):
          edge = geompy.GetSubShape( faceGeomObject, [edges[j]] )
          if geompy.MinDistance( vertex, edge ) == 0 : relatedEdges.append( geompy.GetSubShapeID(faceGeomObject, edge) )
       SE[vertices[i]] = relatedEdges

    #Creating and filling the list of ordered vertices
    continuousVertices = []
    """
    The first vertex added to the list of continuous vertices will be the first key of SE then SE is parsed and :
    1 - for each vertex, getting the two corresponding edges
    2 - taking one of these two corresponding edges and getting its correspondant vertices in ES
    3 - one of these two corersponding vertices is necessarily the initial vertex and will be eliminated, the second one will be added to the list of continuous vertices and so on
    """
    allVertices = list(SE.keys())
    initialVertex = allVertices[0]
    initialEdge = SE[initialVertex][0]
    #Initialization
    currentVertex = initialVertex
    currentEdge = initialEdge
    #Parsing all the edges
    nextVertex = None
    while nextVertex != initialVertex:
       continuousVertices.append( currentVertex )

       #Getting the next edge
       relatedEdges = SE[currentVertex]
       edge_1 = relatedEdges[0]
       edge_2 = relatedEdges[1]
       if edge_1 != currentEdge: 
         nextEdge = edge_1
       else:
         nextEdge = edge_2

       #Getting the next vertex
       relatedVertices = ES[nextEdge]
       vertex_1 = relatedVertices[0]
       vertex_2 = relatedVertices[1]
       if vertex_1 != currentVertex : 
         nextVertex = vertex_1
       else:
         nextVertex = vertex_2

       currentVertex = nextVertex
       currentEdge = nextEdge

    #Transforming the elements of the resulting list from IDs to geom_objects
    for i in range( len(continuousVertices) ):
       continuousVertices[i] = geompy.GetSubShape( faceGeomObject, [continuousVertices[i]] )
    return continuousVertices

def getConvexOrConcaveEdges( shape ):
    geompy = XSS.getGeompy()
    epsilon = 10e-6
    #Global lists of convex/concave edges related to the shape(solid)
    globalConvexList = []
    globalConcaveList = []
    faces = geompy.SubShapeAll( shape, geompy.ShapeType['FACE'] )
    for i in range( len(faces) ):
       #Local lists of convex/concave edges related to each face
       convexList = []
       concaveList = []
       convexVertices = []
       concaveVertices = []
       face = faces[i]
       vertices = getContinuousListOfVertices( face )
       vertices.append( vertices[0] )
       vertices.append( vertices[1] )
       edges = geompy.SubShapeAll( face, geompy.ShapeType['EDGE'] )

       #Dico of correspondances vertices-->edges
       SE = {}
       for i in range( len(vertices) ):
          vertex = vertices[i]
          relatedEdges = []
          for j in range( len(edges) ) :
             edge = edges[j]
             if geompy.MinDistance( vertex, edge ) == 0 : 
                relatedEdges.append( geompy.GetSubShapeID(face, edge) )
          SE[geompy.GetSubShapeID(face, vertex)] = relatedEdges

       parser = 0
       while parser < len(vertices) -2:
          s0 = vertices[parser]
          s1 = vertices[parser+1]
          s2 = vertices[parser+2]
          edge = geompy.MakeEdge( s0, s2 )
          center = geompy.MakeCDG( edge )
          distanceToFace = geompy.MinDistance( center, face )
          if distanceToFace < epsilon : 
            convexVertices.append( s1 )
          else:
            concaveVertices.append( s1 )
          parser += 1

       for i in range( len(convexVertices) ):
          vertex = convexVertices[i]
          id = geompy.GetSubShapeID( face, vertex )
          relatedEdgesIDs = SE[id]
          convexList.append( relatedEdgesIDs[0] )
          convexList.append( relatedEdgesIDs[1] )

       for i in range( len(concaveVertices) ):
          vertex = concaveVertices[i]
          id = geompy.GetSubShapeID( face, vertex )
          relatedEdgesIDs = SE[id]
          concaveList.append( relatedEdgesIDs[0] )
          concaveList.append( relatedEdgesIDs[1] )

       convexList = list( set(convexList) - set(concaveList) )
       #Removing duplicated elements
       convexList = list( set(convexList) )
       concaveList = list( set(concaveList) )
       #Transforming the elements from IDs to geom_objects
       for i in range( len(convexList) ):
          convexList[i] = geompy.GetSubShape( face, [convexList[i]] )

       for i in range( len(concaveList) ) :
          concaveList[i] = geompy.GetSubShape( face, [concaveList[i]] )

       globalConvexList += convexList
       globalConcaveList += concaveList

    #Transforming the elements of the global lists from geom_objects to IDs
    for i in range( len(globalConvexList) ):
      globalConvexList[i] = geompy.GetSubShapeID( shape, globalConvexList[i] )
    for i in range( len(globalConcaveList) ):
      globalConcaveList[i] = geompy.GetSubShapeID( shape, globalConcaveList[i] )

    globalConvexList = list( set(globalConvexList) - set(globalConcaveList) )
    #Removing duplicated elements
    globalConvexList = list( set(globalConvexList) )
    globalConcaveList = list( set(globalConcaveList) )
    #Transforming the elements from IDs to geom_objects
    for i in range( len(globalConvexList) ):
       globalConvexList[i] = geompy.GetSubShape( shape, [globalConvexList[i]] )

    for i in range( len(globalConcaveList) ):
       globalConcaveList[i] = geompy.GetSubShape( shape, [globalConcaveList[i]] )

    return [globalConvexList,globalConcaveList]
