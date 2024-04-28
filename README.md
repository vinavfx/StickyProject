<b>Original Source:</b> https://www.nukepedia.com/plugins/3d/stickyproject
by Ivan Busquets

## StickyProject:
StickyProject is a modified version of UVProject that allows to 'freeze' the projection at a given frame, and let the texture 'stick' to animated geometry from there on.

## Workflow:
Using StickyProject is very similar to using UVProject. You connect it to a camera and a geometry/scene, and it modifies the UVs of that geometry based on the frustum of that camera. A vertex that falls exactly on the lower-left corner of the camera view will get a UV coordinate of 0,0, and a vertex that falls on the top-right corner gets a coordinate of 1,1.

## Compile Source:
```sh
# requires gcc 6.3 on centos 7
cd ~/.nuke/StickyProject
gcc -shared -fPIC -I/opt/Nuke<version>/include -o build/StickyProject.so StickyProject.cpp
```

## Add to menu.py
```python

```


