/* Fungsi pencarian SVD dengan Algoritma Golub-Reinsch 
 * Referensi: Singular value decomposition and least squares solutions
              https://doi.org/10.1007/BF02163027
              Tambahan: http://statistics.uchicago.edu/~lekheng/courses/324/chan.pdf */
function svdoptimized(m1) {
    const withu = true;
    const withv = true;
    let eps = 1e-15;
    let tol = Number.MIN_VALUE / eps;
    var i, j, a, nn, s, h, l, M, d, p, b, u, w;
    let m = m1.length;
    let n = m1[0].length;
    let istransposed = false;
    // Melakukan transpose jika n > m
    if(n > m){
        istransposed = true;
        m1 = math.transpose(math.matrix(m1))._data;
        m = m1.length;
        n = m1[0].length;
    }
    
    var v = [];
    var U = [];
    var VT = [];
    var b = M = 0;
    for (i = 0; i < m; i++)
        U[i] = createzeroarray(m);
    for (i = 0; i < n; i++)
            VT[i] = createzeroarray(n);
    var S = createzeroarray(n);
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            U[i][j] = m1[i][j];

    for (i = 0; i < n; i++) {
        for (v[i] = M, p = 0, nn = i + 1, j = i; j < m; j++) p += Math.pow(U[j][i], 2);
        if (p < tol) M = 0;
        else
            for (d = (l = U[i][i]) * (M = l < 0 ? Math.sqrt(p) : -Math.sqrt(p)) - p, U[i][i] = l - M, j = nn; j < n; j++) {
                for (p = 0, a = i; a < m; a++) p += U[a][i] * U[a][j];
                for (l = p / d, a = i; a < m; a++) U[a][j] = U[a][j] + l * U[a][i]
            }
        for (S[i] = M, p = 0, j = nn; j < n; j++) p += Math.pow(U[i][j], 2);
        if (p < tol) M = 0;
        else {
            for (d = (l = U[i][i + 1]) * (M = l < 0 ? Math.sqrt(p) : -Math.sqrt(p)) - p, U[i][i + 1] = l - M, j = nn; j < n; j++) v[j] = U[i][j] / d;
            for (j = nn; j < m; j++) {
                for (p = 0, a = nn; a < n; a++) p += U[j][a] * U[i][a];
                for (a = nn; a < n; a++) U[j][a] = U[j][a] + p * v[a]
            }
        }
        b < (u = Math.abs(S[i]) + Math.abs(v[i])) && (b = u)
    }

    // Accumulation of right-hand transformations
    if (withv)
        for (i = n - 1; 0 <= i; i--) {
            if (0 !== M) {
                for (d = U[i][i + 1] * M, j = nn; j < n; j++) VT[j][i] = U[i][j] / d;
                for (j = nn; j < n; j++) {
                    for (p = 0, a = nn; a < n; a++) p += U[i][a] * VT[a][j];
                    for (a = nn; a < n; a++) VT[a][j] = VT[a][j] + p * VT[a][i]
                }
            }
            for (j = nn; j < n; j++) VT[i][j] = 0, VT[j][i] = 0;
            VT[i][i] = 1, M = v[i], nn = i
        }
    
    // Accumulation of left-hand transformations
    if (withu) {
        for (i = n; i < m; i++) {
            for (j = n; j < m; j++) U[i][j] = 0;
            U[i][i] = 1
        }
        for (i = n - 1; i >= 0; i--) {
            nn = i + 1;
            M = S[i];
            for (j = nn; j < m; j++)
                U[i][j] = 0;
            if (M !== 0) {
                d = U[i][i] * M;
                for (j = nn; j < m; j++) {
                    for (p = 0, a = nn; a < m; a++) p += U[a][i] * U[a][j];
                    for (l = p / d, a = i; a < m; a++) U[a][j] = U[a][j] + l * U[a][i]
                }
                for (j = i; j < m; j++) U[j][i] = U[j][i] / M
            } else
                for (j = i; j < m; j++) U[j][i] = 0;
            U[i][i] = U[i][i] + 1
        }
    }

    // Diagonalization of the bidiagonal form
    eps *= b;
    for (a = n - 1; a >= 0; a--){
        for (var k = 0; k < 50; k++) {
            let isconvergence = false;
            for (nn = a; nn >= 0; nn--) {
                if (Math.abs(v[nn]) <= eps) {
                    isconvergence = true;
                    break;
                }
                if (Math.abs(S[nn - 1]) <= eps)
                    break;
            }
            // cancellation
            if (!isconvergence)
                for (h = 0, s = nn - (p = 1), i = nn; i < a + 1 && (l = p * v[i], v[i] = h * v[i], !(Math.abs(l) <= eps)); i++){
                    if (M = S[i], S[i] = Math.sqrt(l * l + M * M), h = M / (d = S[i]), p = -l / d, withu){
                        for (j = 0; j < m; j++){
                            u = U[j][s];
                            w = U[j][i];
                            U[j][s] = u * h + w * p;
                            U[j][i] = -u * p + w * h;
                    }
                }
            }
            // test f convergence
            w = S[a];
            if (nn === a) {
                // convergence
                if (w < 0){
                    S[a] = -w;
                    if(withv)
                        for (j = 0; j < n; j++) VT[j][a] = -VT[j][a];
                }
                break;
            }
            b = S[nn];
            l = (((u = S[a - 1]) - w) * (u + w) + ((M = v[a - 1]) - (d = v[a])) * (M + d)) / (2 * d * u);
            M = Math.sqrt(l * l + 1);
            l = ((b - w) * (b + w) + d * (u / (l < 0 ? l - M : l + M) - d)) / b;
            for (i = nn + (p = h = 1); i < a + 1; i++) {
                if (M = v[i], u = S[i], d = p * M, M *= h, w = Math.sqrt(l * l + d * d), l = b * (h = l / (v[i - 1] = w)) + M * (p = d / w), M = -b * p + M * h, d = u * p, u *= h, withv)
                    for (j = 0; j < n; j++) b = VT[j][i - 1], w = VT[j][i], VT[j][i - 1] = b * h + w * p, VT[j][i] = -b * p + w * h;
                if (w = Math.sqrt(l * l + d * d), l = (h = l / (S[i - 1] = w)) * M + (p = d / w) * u, b = -p * M + h * u, withu)
                    for (j = 0; j < m; j++) u = U[j][i - 1], w = U[j][i], U[j][i - 1] = u * h + w * p, U[j][i] = -u * p + w * h
            }
            v[nn] = 0;
            v[a] = l;
            S[a] = b;
        }
    }

    for (i = 0; i < n; i++){
        if(S[i] < eps)
            S[i] = 0;
    }

    if(istransposed){   // Membalikkan u dan v jika n > m
        return {u: VT, q: S, v: U};
    } else{
        return {u: U, q: S, v: VT};
    }
}

function createzeroarray(n){
    var arr = new Array(n);
    for(let i = 0; i < n; i++){
        arr[i] = 0;
    }
    return arr;
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

function compress(m1){
    let mtrx = [];
    for (let i = 0; i < m1.rows; i++) {
        mtrx.push(Array.from(m1.data.slice(i*m1.cols, ((i+1)*m1.cols))));
    }

    let svd = svdoptimized(mtrx);

    let m = mtrx.length; // baris
    let n = mtrx[0].length; // kolom

    const mtrx_zero = new Array(m).fill(0).map(() => new Array(n).fill(0));

    mtrx_zero[0][0] = svd.q[0];
    for (var i = 0; i < Math.min(m, n); i++){
        for (var j = 0; j < Math.min(m, n); j++){
            if (i = j){
                mtrx_zero[i][j] = svd.q[i];
            }
        }
    }

    size_compress = 100; // ini buat tingkat kompresinya

    const mtrx_u = new Array(m).fill(0).map(() => new Array(size_compress).fill(0));

    for (var i = 0; i < m; i++){
        for (var j = 0; j < size_compress; j++){
            mtrx_u[i][j] = svd.u[i][j];
        }
    }

    const mtrx_q = new Array(size_compress).fill(0).map(() => new Array(size_compress).fill(0));

    for (var i = 0; i < size_compress; i++){
        for (var j = 0; j < size_compress; j++){
            mtrx_q[i][j] = mtrx_zero[i][j];
        }
    }

    const mtrx_v = new Array(size_compress).fill(0).map(() => new Array(n).fill(0));

    for (var i = 0; i < size_compress; i++){
        for (var j = 0; j < n; j++){
            mtrx_v[i][j] = svd.v[i][j];
        }
    }

    let m2 = [];
    m2 = math.multiply(mtrx_u, math.multiply(mtrx_q, mtrx_v));

    return m2;
}
