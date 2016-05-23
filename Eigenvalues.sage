
def gersgorin_circle (m) :                                                                     # gersgorin circle
    c = circle ((0,0),0)
    for i in range (0, m.nrows()) :
        r = 0
        for k in range (0, m.ncols()) :
            if i == k : continue
            r = r + abs(m[i,k])
        c = c + circle ((m[i,i], 0), r)

    show (c)

G = matrix (QQ, [[4,1,1], [0,2,1], [-2,0,9]])

def dot_ (a, b) :                                                                             # dot product
    if a.nrows() != b.nrows() :
        return "Error"
    ans = 0
    for x,y in zip (a,b) :
        ans = ans + x * y
    return ans

def magnitude (v) :                                                                          # magnitude of a vector
    ans = 0
    for x in v :
        ans = ans + x*x
    return sqrt (ans)

def gram_schmidt (m) :                                                                             # gram schmidt process

    nc = m.ncols()
    A = m[:,0]/magnitude(m[:,0])     # The first column of 'm' is the first vector in our orthogonal set A

    for i in range (1, nc) :
        v = m[:,i]  
        for k in range (0,i) :
            d = dot_(A[:,k], m[:,i]) / dot_(A[:,k], A[:,k])
            v = v - d * A[:,k]
        v = v / magnitude (v)                 ## HAVING PROBLEMS WITH TYPES.
        A = A.augment (v)

    return A

A = matrix (QQ, [[1,1,1], [0,1,1], [0,0,1]])
B = matrix (RDF, [[0,-5,1], [4,-1,-1], [2,2,2]])

def is_orthogonal (m) :                                                                     # determine if a matrix is orthogonal
    # any invertible matrix m with m^-1 = m^t is orthogonal
    # Also, if m*m^t == I, then m^t = m^-1.

    if m.nrows() != m.ncols() :
        return False

    if m * m.transpose() == identity_matrix (m.nrows()) :
        return True
    return False

def power_method (A, x, TOL, N) :                                             # matrix A, intial guess vector x, Tolerance TOL, Number of iterations N
    k = 1
    p = 0
    for cur in range (0,len(x)) :
        if abs(x[cur]) > abs(x[p]) :
            p = cur
    x = x / x[p]

    while k <= N :
        y = A * x
        u = y[p]

        p = 0
        for cur in range (0, len(y)) :
            if abs (y[cur]) > abs (y[p]) :
                p = cur

        if y[p] == 0 : return False, "A has eigenvalue 0. Select a new vector x and restart"

        ERR = 0
        for x1, y1 in zip(x,y) :
            ERR = max (ERR, x1 - y1 / y[p])
        x = y / y[p]

        if ERR < TOL :
            return u, x

        k = k + 1

    return False, "The maximum number of iterations exceeded"

A_ = matrix (QQ, [[-2,-3], [6,7]])
x_A = vector (QQ, [1,1])


def householder (A, k) :

    n = A.nrows()

    q = 0
    for j in range (k+1, n) :
        q = q + A[j,k] ** 2
    q = q ** 0.5

    alpha = -q
    if A[k+1,k] != 0:
        alpha = -q * A[k+1,k] / abs (A[k+1,k])

    r = (0.5 * alpha**2 - 0.5*alpha*A[k+1,k]) ** 0.5

    w = vector (RR, {(n-1): 0})
    w[k+1] = (A[k+1,k] - alpha) / (2*r)
    for j in range (k+2, n) :
        w[j] = A[j,k] / (2*r)

    ww_t = matrix (RR, n,n, 0)
    for i in range (0,n) :
        for k in range (0,n) :
            ww_t[i,k] = w[i] * w[k]

    P = identity_matrix (n) - 2 * ww_t
    A = P * A * P

    k = k + 1

    if k >= n-1 :
        return A
    householder(A, k)

#AA = matrix (RR, [[4,1,-2,2], [1,2,0,1], [-2,0,3,-2], [2,1,-2,-1]])


def QR (A) :
    
    Q = gram_schmidt(A)
    R = Q.transpose() * A
    
    return Q, R


M = matrix (QQ, [[12,-51,4], [6, 167, -68], [-4, 24, -41]])

def least_squares (A, b) :     # Finds approximate solution to Ax = b
    
    Q, R = QR (A)
    xx = R.inverse() * Q.transpose() * b
    return xx

AA_ = matrix (RDF, [[1,3,5], [1,1,0], [1,1,2], [1,3,3]])
b = matrix (RDF, [[3],[5],[7],[-3]])
