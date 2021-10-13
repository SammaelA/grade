bl_info = {
    "name": "TreeGen",
    "category": "Object",
    "description": "Generate high quality tree models",
    "author": "Charlie Hewitt and Luke Pflibsen-Jones",
    "version": (0, 0, 3),
    "wiki_url": "https://github.com/friggog/tree-gen/wiki",
    "tracker_url": "https://github.com/friggog/tree-gen/issues/new/choose",
    'blender': (2, 80, 0)
}


from . import gui


def register():
    print("reg")


def unregister():
    print("unreg")
    # Reversing order is best-practice


if __name__ == "__main__":
    register()
