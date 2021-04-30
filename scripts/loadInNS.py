#loadInNS.py
def load(file):
    exec(open(file).read())
    return locals()
