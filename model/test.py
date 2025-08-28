from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag

## https://github.com/tanghaibao/goatools/blob/main/notebooks/parents_and_ancestors.ipynb

GO_ID = 'GO:0019222' 

godag = GODag('data/GO/go-basic.obo',
              optional_attrs={'relationship'})

#print(godag)

gosubdag_r0 = GoSubDag([GO_ID], godag, prt=None)
print('{GO} ancestors: {P}'.format(
    GO=GO_ID,
    P=gosubdag_r0.rcntobj.go2ancestors[GO_ID]))

def get_ancestors(goid, godag):
    gosubdag_r0 = GoSubDag([GO_ID], godag, prt=None)