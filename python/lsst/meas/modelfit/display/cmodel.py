class CModelDisplay(object):

    def __init__(self, schema, refSchema=None, baseName="modelfit_CModel"):
        self.schema = schema
        self.refSchema = refSchema
        self.ellipseKeys = {}
        self.fluxKeys = {}
        for component in ("initial", "exp", "dev"):
            prefix = "_".join((baseName, component))
            self.fluxKeys[component] = schema.find("{}_fluxInner".format(prefix)).key
            if self.refSchema is None:
                self.ellipseKeys = {schema.find("_".join([baseName] + ))}



def getComponentShapelets(record, name, reference=None):
