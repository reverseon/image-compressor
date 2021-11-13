/* Fungsi pencarian SVD dengan Algoritma Golub-Reinsch 
 * Referensi: Singular value decomposition and least squares solutions
              https://doi.org/10.1007/BF02163027
              Tambahan: http://statistics.uchicago.edu/~lekheng/courses/324/chan.pdf */
function svdgolub(m1) {
    withu = true;
    withv = true;
    let m = m1.length;
    let n = m1[0].length;
    let eps = 1e-15;
    let tol = Number.MIN_VALUE / eps;

    let U;
    let istransposed = false;
    if(n > m){
        istransposed = true;
        U = matrixtranspose(m1);
        m = U.length;
        n = U[0].length;
    } else{
        U = matrixclone(m1);
    }
    let Q = createzeroarray(n);
    let V = createzeromatrix(n, n);
    let e = createzeroarray(n);

    let i, j, k, l, l1;
    let c, f, g, h, s, x, y, z;

    // Reduksi Householder ke bentuk bidiagonal
    g = 0;
    x = 0;
    for (i = 0; i < n; i++) {
        e[i] = g;
        s = 0;
        l = i + 1;
        for (j = i; j < m; j++) {
            s = s + U[j][i] ** 2;
        }
        if (s < tol) {
            g = 0;
        } else {
            f = U[i][i];
            g = (f < 0) ? (s ** 0.5) : -(s ** 0.5);
            h = f * g - s;
            U[i][i] = f - g;

            for (j = l; j < n; j++) {
                s = 0;
                for (k = i; k < m; k++) {
                    s = s + U[k][i] * U[k][j];
                }
                f = s / h;
                for (k = i; k < m; k++) {
                    U[k][j] = U[k][j] + f * U[k][i];
                }
            }
        }
        Q[i] = g;
        s = 0;
        for (j = l; j < n; j++) {
            s = s + U[i][j] ** 2;
        }
        if (s < tol)
            g = 0;
        else {
            f = U[i][i + 1];
            g = (f < 0) ? (s ** 0.5) : -(s ** 0.5);
            h = f * g - s;
            U[i][i + 1] = f - g;
            for (j = l; j < n; j++) {
                e[j] = U[i][j] / h;
            }
            for (j = l; j < m; j++) {
                s = 0;
                for (k = l; k < n; k++) {
                    s = s + U[j][k] * U[i][k];
                }
                for (k = l; k < n; k++) {
                    U[j][k] = U[j][k] + s * e[k];
                }
            }
        }
        y = Math.abs(Q[i]) + Math.abs(e[i]);
        if (y > x) x = y;
    }

    // Accumulation of right-hand transformations
    if (withv) {
        for (i = n - 1; i >= 0; i--) {
            if (g != 0) {
                h = U[i][i + 1] * g;
                for (j = l; j < n; j++) V[j][i] = U[i][j] / h;
                for (j = l; j < n; j++) {
                    s = 0;
                    for (k = l; k < n; k++) s = s + U[i][k] * V[k][j];
                    for (k = l; k < n; k++) V[k][j] = V[k][j] + s * V[k][i];
                }
            }
            for (j = l; j < n; j++) V[i][j] = V[j][i] = 0;
            V[i][i] = 1;
            g = e[i];
            l = i;
        }
    }

    // Accumulation of left-hand transformations
    if (withu) {
        for (i = n; i < m; i++) {
            for (j = n; j < m; j++) U[i][j] = 0;
            U[i][i] = 1;
        }
        for (i = n - 1; i >= 0; i--) {
            l = i + 1;
            g = Q[i];
            for (j = l; j < m; j++) U[i][j] = 0;
            if (g != 0) {
                h = U[i][i] * g;
                for (j = l; j < m; j++) {
                    s = 0;
                    for (k = l; k < m; k++) s += U[k][i] * U[k][j];
                    f = s / h;
                    for (k = i; k < m; k++) U[k][j] += f * U[k][i];
                }
                for (j = i; j < m; j++) U[j][i] = U[j][i] / g;
            } else {
                for (j = i; j < m; j++) U[j][i] = 0;
            }
            U[i][i] = U[i][i] + 1;
        }
    }
    
    // Diagonalization of the bidiagonal form
    eps = eps * x;
    k = n - 1;
    let counter = 0;
    //for (k = n - 1; k >= 0; k--) {
    while(k >= 0){
        // test f splitting:
        let cancel = true;
        for (l = k; l >= 0; l--) {
            if (Math.abs(e[l]) <= eps){
                // goto test f convergence (do nothing)
                cancel = false;
                break;
            } 
            if (Math.abs(Q[l - 1]) <= eps){
                cancel = true; // goto cancellation
                break;
            }
        }

        if (cancel) {
            // cancellation:
            c = 0;
            s = 1;
            l1 = l - 1;
            for (i = l; i < k; i++) {
                f = s * e[i];
                e[i] = c * e[i];
                if (Math.abs(f) <= eps) {
                    break;
                } // goto test f convergence
                g = Q[i];
                h = Q[i] = Math.sqrt(f * f + g * g);
                c = g / h;
                s = -f / h;
                if (withu) {
                    for (j = 0; j < m; j++) {
                        y = U[j][l1];
                        z = U[j][i];
                        U[j][l1] = y * c + z * s;
                        U[j][i] = -y * s + z * c;
                    }
                }
            }
        }
        
        // test f convergence:
        z = Q[k];
        if (l == k) {
            //goto convergence
            // convergence:
            if (z < 0) {
                Q[k] = -z;
                if (withv) {
                    for (j = 0; j < n; j++) V[j][k] = -V[j][k];
                }
            }
            counter = 0;
            k--;
        } else {
            // Shift from bottom 2x2 minor
            x = Q[l];
            y = Q[k - 1];
            g = e[k - 1];
            h = e[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
            g = Math.sqrt(f * f + 1);
            f = ((x - z) * (x + z) + h * (y / (f < 0 ? f - g : f + g) - h)) / x;

            // Next QR transformation
            c = 1;
            s = 1;
            for (i = l + 1; i < k + 1; i++) {
                g = e[i];
                y = Q[i];
                h = s * g;
                g = c * g;
                e[i - 1] = z = Math.sqrt(f * f + h * h);
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = -x * s + g * c;
                h = y * s;
                y = y * c;
                if (withv) {
                    for (j = 0; j < n; j++) {
                        x = V[j][i - 1];
                        z = V[j][i];
                        V[j][i - 1] = x * c + z * s;
                        V[j][i] = -x * s + z * c;
                    }
                }
                Q[i-1] = z = Math.sqrt(f*f + h*h);
                c = f/z;
                s = h/z;
                f = c*g + s*y;
                x = -s*g + c*y;
                if (withu) {
                    for (j = 0; j < m; j++) {
                        y = U[j][i - 1];
                        z = U[j][i];
                        U[j][i - 1] = y * c + z * s;
                        U[j][i] = -y * s + z * c;
                    }
                }
            }
            e[l] = 0;
            e[k] = f;
            Q[k] = x;
            // GOTO test f splitting (loop)
            counter++;
            if(counter == 5){
                k--;
                counter = 0;
            }
        }
    }
    for (i = 0; i < n; i++){
        if(Q[i] < eps)
            Q[i] = 0;
    }
    
    if(istransposed){
        return {u: V, q: Q, v: U};
    } else{
        return {u: U, q: Q, v: V};
    }
}

function matrixclone(m1) {
    let m2 = [];
    m = m1.length;
    n = m1[0].length;
    for (let i = 0; i < m; i++) {
        let temp = []
        for (let j = 0; j < n; j++) {
            temp.push(m1[i][j]);
        }
        m2.push(temp);
    }
    return m2;
}

function matrixtranspose(m1){
    if(m1.length == 0 || m1[0].length == 0)
        return [];
    let m = m1.length;
    let n = m1[0].length;
    let m2 = []
    for(let i=0; i<n; i++){
        let temp = []
        for(let j=0; j<m; j++){
            temp.push(m1[j][i]);
        }
        m2.push(temp);
    }
    return m2;
}

function createzeroarray(n){
    var arr = new Array(n);
    for(let i = 0; i < n; i++){
        arr[i] = 0;
    }
    return arr;
}

function createzeromatrix(m, n) {
    let m1 = [];
    for (let i = 0; i < m; i++) {
        let temp = [];
        for (let j = 0; j < n; j++) {
            temp.push(0);
        }
        m1.push(temp);
    }
    return m1;
}