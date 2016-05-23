
def rref (m) : # row reduced row echelon form. input a matrix

    tolerance = 0.00000000001     # Round-off errors.

    p = 0
    nrow = m.nrows()
    ncol = m.ncols()

    for r in range (0, nrow) :

        while True :     # go through columns until we find a pivot column

            r1 = r
            while r1 < nrow and abs(m[r1][p]) <= tolerance :         # try the current column to find nonzero value.
                r1 = r1 + 1                                     # (or a value greater than the "tolerance")

            if r1 >= nrow :             # The current row has no non-zero entry
                p = p + 1               # so look at the next column.
            else :
                break

            if p >= ncol :         # we have reached past the last column, so we're done.
                return m

        m.swap_rows(r,r1)
        pivot_ = m[r,p]
        for c in range (p,ncol) :
            m[r,c] = m[r,c] / pivot_

        for j in range (0, nrow) :
            if j != r :
                m[j] = m[j] - m[j,p] * m[r]


        p = p + 1              # move on to next column
        if p >= ncol :         # if we have inspected each column, we are done.
            break

    return m


# TEST MATRICES. -- it works

A = matrix (QQ, [[0,3,-6,6,4,-5],[3,-7,8,-5,8,9],[3,-9,12,-9,6,15]])
B = matrix (QQ, [[1,-2],[-2,7],[3,-9]])
C = matrix (QQ, [[1,6,2,-4],[-3,2,-2,-8],[4,-1,3,9]])
D = matrix (QQ, [[1,14,2],[23,2,66],[1,2910,-14]])
E = matrix (QQ, [[-3,6,-1,1,-7],[1,-2,2,3,-1],[2,-4,5,8,-4]])   # Exemplifies round off error (if you set tolerance to 0, we get a wrong answer)

def solve_system (m) :     # system of linear equations with variables x_1, x_2, ... x_n

    nrow = m.nrows()
    ncol = m.ncols()

    m = rref (m)

    # find any free variables. -- There would be a row of zeros

    free = []
    bound = []

    for r in range (0,nrow) :             # first, we determine the bound variables
        for c in range (0,ncol) :
            if c == ncol-1 and m[r,c] != 0 :     # last column is a pivot column, so the system is inconsistent
                return "No Solution"
            if m[r,c] != 0 :
                bound.append(c)
                break

    for c in range (0,ncol-1) :          # if a variable's not bound, it must be freee.
        if c not in bound :
            free.append (c)

     # print out answer
    i = 0
    for b in bound :

        print "x_%d = " %b,
        if m[i, ncol-1] != 0 :
            print "%d" %m[i,ncol-1],

        for f in free :
            if m[i,f] > 0 :
                print "- %dx_%d" %(m[i,f] ,f),
            elif m[i,f] < 0 :
                 print "+ %dx_%d" %(abs(m[i,f]) ,f),

        i = i + 1
        print

    if len(free) > 0 :
        print "Where variables",
        for f in free :
            print "x_%d," %f,
        print "are free."
        print


AA = matrix (QQ, [[1,-1,2,-1,-8], [2,-2,3,-3,-20], [1,1,1,0,-2], [1,-1,4,3,4]])    # no free variables
BB = matrix (QQ, [[1,-3,-5,0], [0,1,1,3]])                                # one free variable
CC = matrix (QQ, [[1,-2,-1,3,0], [-2,4,5,-5,3],[3,-6,-6,8,3]])        # inconsisten system
DD = matrix (QQ, [[1,6,2,-5,-2,-4],[0,0,2,-8,-1,3], [0,0,0,0,1,7]])    # two free variables.

def inverse (m) :

    nrow = m.nrows()
    ncol = m.ncols()
    if nrow != ncol :
        return "Noninvertible Matrix"

    I = identity_matrix(nrow)
    m = m.augment (I)
    rref (m)

    left = m.submatrix(0,0,nrow,ncol)
    right = m.submatrix(0,ncol, nrow, ncol)

    if left != I :
        return "Nonsingular Matrix"
    return right

A_ = matrix (QQ, [[1,2,-1], [2,1,0], [-1,1,2]])      # invertible matrix
B_ = matrix (QQ, [[1,-3,-5,0], [0,1,1,3]])           # Noninvertible
C_ = matrix (QQ, [[1,14,2],[23,2,0],[1,2910,-14]])   # nonsigular

def LU (m) :
    tolerance = 0.00000000001

    p = 0
    nrow = m.nrows()
    ncol = m.ncols()

    L = identity_matrix (QQ, nrow)    # lower identiy matrix
    Per = identity_matrix (nrow)    # permutation matrix

    # this is the same as in rref
    for r in range (0, nrow) :

        while True :
            r1 = r
            while r1 < nrow and abs(m[r1][p]) <= tolerance :
                r1 = r1 + 1
            if r1 >= nrow :
                p = p + 1
            else :
                break
            if p >= ncol :
                return m, L, Per

        m.swap_rows(r,r1)
        Per.swap_rows(r,r1)

        pivot_ = m[r,p]
        for rr in range (r,nrow) :
            L[rr,r] = m[rr,p] / pivot_

        for j in range (r, nrow) :
            if j != r :
                temp = m[j,p]
                m[j] = m[j] - (temp / m[r,p]) * m[r]

        p = p + 1              # move on to next column
        if p >= ncol :         # if we have inspected each column, we are done.
            break
    return Per, L, m

AA_ = matrix (QQ, [[1,1,0,3], [2,1,-1,1], [3,-1,-1,2], [-1,2,3,-1]])
BB_ = matrix (QQ, [[2,4,-1,5,-2], [-4,-5,3,-8,1], [2,-5,-4,1,8], [-6,0,7,-3,1]])

def determinant (m) :

    nrow = m.nrows()
    ncol = m.ncols()

    if nrow != ncol :
        return "This is not an invertible matrix. Therefore, it does not have a determinant."

    P,L,U = LU (m)

    d = 1
    for i in range (0,nrow) :
        d = d * L[i,i]
        d = d * U[i,i]

    # find how many row swaps it takes to get from P to identity. If that number is even, det(P) = 1, else det(P) = -1
    swaps = 0
    for i in range (0,nrow) :
        if P[i,i] == 1 : continue
        for k in range (i+1,nrow) :
            if P[i,k] == 1 :
                P[k] = P[i]
                swaps = swaps + 1

    if swaps % 2 == 0 :
        return d
    return d * (-1)


A_d = matrix (QQ, [[2,1,-1,1], [1,1,0,3], [-1,2,3,-1], [3,-1,-1,2]])
B_d = matrix (QQ, [[2,-1,3,0], [4,-2,7,0], [-3,-4,1,5], [6,-6,8,0]])









