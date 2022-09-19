import drjit as dr
import mitsuba as mi
from matplotlib import pyplot as plt
import csv
from mitsuba.scalar_rgb import Transform4f as T
import struct

def init(base_path):
  print("base_path = ", base_path)
  mi.set_variant('cuda_ad_rgb')
  integrator = {
      'type': 'direct_reparam',
  }
  scene = mi.load_dict({
      'type': 'scene',
      'integrator': integrator,
      'sensor':  {
          'type': 'perspective',
          'to_world': T.look_at(
                          origin=(0, 0, 2),
                          target=(0, 0, 0),
                          up=(0, 1, 0)
                      ),
          'fov': 60,
          'film': {
              'type': 'hdrfilm',
              'width': 256,
              'height': 256,
              'rfilter': { 'type': 'gaussian' },
              'sample_border': True
          },
      },
      'wall': {
          'type': 'obj',
          'filename': base_path + 'meshes/rectangle.obj',
          'to_world': T.translate([0, 0, -2]).scale(2.0),
          'face_normals': True,
          'bsdf': {
              'type': 'diffuse',
              'reflectance': { 'type': 'rgb', 'value': (0.5, 0.5, 0.5) },
          }
      },
      'bunny': {
          'type': 'obj',
          'filename': base_path + 'meshes/sphere.obj',
          'to_world': T.scale(0.25),
          'bsdf': {
              'type': 'diffuse',
              'reflectance': { 'type': 'rgb', 'value': (0.3, 0.3, 0.75) },
          },
      },
      'light': {
          'type': 'obj',
          'filename': base_path + 'meshes/sphere.obj',
          'emitter': {
              'type': 'area',
              'radiance': {'type': 'rgb', 'value': [1e3, 1e3, 1e3]}
          },
          'to_world': T.translate([2.5, 2.5, 7.0]).scale(0.25)
      }
  })

  img_ref = mi.render(scene, seed=0, spp=256)

  mi.util.convert_to_bitmap(img_ref)
  params = mi.traverse(scene)
  initial_vertex_positions = dr.unravel(mi.Point3f, params['bunny.vertex_positions'])
  opt = mi.ad.Adam(lr=0.025)
  opt['angle'] = mi.Float(0.0)
  opt['trans'] = mi.Point3f(0, 0.5, 0)
  #opt['bunny.vertex_positions'] = params['bunny.vertex_positions']

  F_transform(params, opt, initial_vertex_positions)

  context = {
    'scene' : scene,
    'img_ref' : img_ref,
    'params' : params, 
    'initial_vertex_positions' : initial_vertex_positions,
    'vertex_positions' : params['bunny.vertex_positions'],
    'vertex_normals' : params['bunny.vertex_normals'],
    'vertex_texcoords' : params['bunny.vertex_texcoords'],
    'opt' : opt
  }
  return context

def F_transform(params, opt, initial_vertex_positions):
    opt['trans'] = dr.clamp(opt['trans'], -10, 10)
    opt['angle'] = dr.clamp(opt['angle'], -1, 1)
    trafo = mi.Transform4f.translate([opt['trans'].x, opt['trans'].y, opt['trans'].z]).rotate([0, 1, 0], opt['angle'] * 100.0)
    tr_positions = trafo @ initial_vertex_positions
    params['bunny.vertex_positions'] = dr.ravel(tr_positions)
    #params['bunny.vertex_positions'] = opt['bunny.vertex_positions']
    params.update()

def F_loss(img, img_ref):
    loss = dr.sum(dr.sqr(img - img_ref)) / len(img)
    return loss

def opt_iter(it, context):
  scene = context['scene']
  params = context['params']
  opt = context['opt']
  initial_vertex_positions = context['initial_vertex_positions']
  img_ref = context['img_ref']

  F_transform(params, opt, initial_vertex_positions) # Scene = F_transform(input)
  img = mi.render(scene, params, seed=it, spp=32) # image = F_render(scene)
  loss = F_loss(img, img_ref) # loss = F_loss(image)
  dr.backward(loss) # calculate F_loss(F_render(F_transform(input))) gradient

  opt.step()
  print(f"Iteration {it:02d}: error={loss[0]:6f}, angle={opt['angle'][0]:.4f}, trans=[{opt['trans'].x[0]:.4f}, {opt['trans'].y[0]:.4f}]", end='\n')
  mi.util.write_bitmap("iter.png", img)
  it = it+1
  return loss[0]

def render(it, context):
  scene = context['scene']
  params = context['params']
  opt = context['opt']
  img_ref = context['img_ref']

  params['bunny.vertex_positions'] = context['vertex_positions']
  params['bunny.vertex_normals'] = context['vertex_normals']
  params['bunny.vertex_texcoords'] = context['vertex_texcoords']
  params.update()

  dr.enable_grad(params['bunny.vertex_positions'])
  dr.enable_grad(params['bunny.vertex_normals'])
  dr.enable_grad(params['bunny.vertex_texcoords'])
  
  img = mi.render(scene, params, seed=it, spp=4) # image = F_render(scene)
  loss = F_loss(img, img_ref) # loss = F_loss(image)
  dr.backward(loss)

  context['vertex_positions_grad'] = dr.grad(params['bunny.vertex_positions'])
  context['vertex_normals_grad'] = dr.grad(params['bunny.vertex_normals'])
  context['vertex_texcoords_grad'] = dr.grad(params['bunny.vertex_texcoords'])

  print(f"Iteration {it:02d}: error={loss[0]:6f}, angle={opt['angle'][0]:.4f}, trans=[{opt['trans'].x[0]:.4f}, {opt['trans'].y[0]:.4f}]", end='\n')
  mi.util.write_bitmap("iter.png", img)
  it = it+1
  return loss[0]

def get_params(context, key):
  lst = list(context[key])
  #print("Python: sent ",len(lst)," floats to C++")
  #for i in range(10):
  #  print(lst[i])
  return struct.pack('%sf' % len(lst), *lst)

def set_params(context, key, bytes, n):
  tup = struct.unpack('%sf' % n, bytes)
  context[key] = tup
  #print("Python: receive ",len(tup)," floats from C++")
  #for i in range(10):
  #  print(tup[i])

def pow2(val):
  return val*2

def standalone_test():
  base_path = '../resources/mitsuba_data/'
  ctx = init(base_path)
  for i in range(50):
    opt_iter(i, ctx)


#standalone_test()