# Import the library using the alias "mi"
import mitsuba as mi
# Set the variant of the renderer
mi.set_variant('scalar_rgb')
# Load a scene
scene = mi.load_dict(mi.cornell_box())
# Render the scene
img = mi.render(scene)
mi.util.write_bitmap("cbox.png", img)