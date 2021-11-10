class UserException {
    constructor(clasification, msg) {
        this.msg = msg;
        this.clsf = clasification;
    }
}

class EigenPair {
    constructor(egval) {
        this.egvec = math.matrix([[]]) // EIGENVECTOR PER ROW, 1 ROW = 1 EIGENVECTOR
        this.egval = egval
    }
}

class SVDM {
    constructor() {
        this.u = math.matrix([[]])
        this.sigma = math.matrix([[]])
        this.vt = math.matrix([[]])
    }
}
let EPSILON = 0.0000000001

function isEqualF(a, b) {
    if (Math.abs(a-b) < EPSILON) {
        return true;
    } else {
        return false;
    }
}

function isMEq(a, b) {
    for (let i = 0; i < a._size[0]; i++) {
        for (let j = 0; j < a._size[1]; j++) {
            if (!isEqualF(a.get([i,j]), b.get([i,j]))) {
                return false; 
            }
        }
    }
    return true;
}

function qrdecompose(m1) {
    let n = m1.size()[0];
    let I = math.identity(n);

    m1 = m1.clone();
    let result = {
        Q: math.identity(n),
        R: m1
    };

    for (let i = 0; i < n - 1; i++) {
        let columnvec = math.transpose(math.column(result.R, i));
        let normvec = math.divide(columnvec, math.norm(columnvec._data[0]));

        let D = math.norm(normvec._data[0].slice(i));
        if (normvec.get([0, i]) > 0) {
            D *= -1;
        }

        let dk = normvec.get([0, i]);
        let vk = (0.5 * (1 - dk / D)) ** 0.5;
        let p = -1 * D * vk;
        
        let VT = math.matrix(math.zeros([1, n]));
        VT.set([0, i], vk);
        for (let k = i+1; k < n; k++) {
            VT.set([0, k], normvec.get([0, k]) / (2 * p));
        }

        let V = math.transpose(VT);
        let tmp1 = math.multiply(V, VT);
        let tmp2 = math.multiply(tmp1, 2);
        
        let Pk = math.subtract(I, tmp2);
        result.Q = math.multiply(Pk, result.Q);
        result.R = math.multiply(Pk, result.R);
    }
    result.Q = math.transpose(result.Q);
    return result;
}

function eigennoshift(m, step) {
    let a = math.clone(m)
    for (let i = 0; i < step; i++) {
        let qr = math.qr(a);
        let prea = a
        a = math.multiply(qr.R, qr.Q)
        if (isMEq(prea, a)) {
            break;
        }
    }
    return a;
}

function eigenshift(m, step) {
    let a = math.clone(m)
    for (let i = 0; i < step; i++) {
        shift = math.multiply(a.get([a._size[0]-1, a._size[1]-1]),math.identity(math.size(a)))
        let prea = a
        a = math.subtract(a, shift);
        qr = math.qr(a);
        q = qr.Q
        r = qr.R
        a = math.add(math.multiply(r, q),  shift);
        if (i < 3) { // arbitrary
            if (isMEq(prea, a)) {
                a = eigennoshift(m, 5);
                break;
            }
        }
    }
    return a;
}

function reducedREF(ma) {
    let m = math.clone(ma._data);
    let lead = 0;
    let row = m.length
    let col = m[0].length
    for (let r = 0; r < row; r++) {
        if (col <= lead) {
            return m;
        }
        let i = r;
        while (isEqualF(m[i][lead], 0)) {
            i++;
            if (row == i) {
                i = r;
                lead++;
                if (col == lead) {
                    return m;
                }
            }
        }
 
        let tmp = m[i];
        m[i] = m[r];
        m[r] = tmp;
 
        let val = m[r][lead];
        for (let j = 0; j < col; j++) {
            m[r][j] /= val;
        }
 
        for (i = 0; i < row; i++) {
            if (i === r) continue;
            val = m[i][lead];
            for (let j = 0; j < col; j++) {
                m[i][j] -= val * m[r][j];
            }
        }
        lead++;
    }
    return m; // ARRAY, NOT A MATRIX
}

function normalizevect(v) {
    let norm = math.norm(v._data[0], 2)
    v = math.divide(v, norm)
    return v;
}

function beautifyv(v) {
    // V is 2D
    let abslowest = Math.abs(v.get([0, 0]))
    for (let i = 1; i < v._size[1]; i++) {
        if (!isEqualF(Math.abs(v.get([0, i])), 0)) {
            if (abslowest > Math.abs(v.get([0, i]))) {
                abslowest = Math.abs(v.get([0, i]));
            }
        }  
    }
    if (!isEqualF(abslowest, 0)) {
        let mltp = 1/abslowest;
        v = math.multiply(mltp, v) 
    }
    return v
}

function isExistinArray(a, v) {
    for (let i = 0; i < a.length; i++) {
        if (isEqualF(a[i], v)) {
            return true;
        }
    }
    return false;
}

function getEigenPair(m) {
    let aoEP = [] // ARRAY OF EIGENPAIR
    let egvalm = eigenshift(m, 5)
    let egval = math.zeros(egvalm._size[0])
    for (let i = 0; i < egvalm._size[0]; i++) {
        egval.set([i], egvalm.get([i,i]))
    }
    let checkedeigen = []
    for (let i = 0; i < egval._size[0];i++) {
        let egvnow = egval.get([i]);
        if (isExistinArray(checkedeigen, egvnow)) {
            continue;
        } else {
            checkedeigen.push(egvnow)
            let tap = new EigenPair(egvnow);
            let adjusting = math.identity(m._size[0])
            adjusting = math.multiply(egvnow, adjusting);
            let tbs = math.subtract(m, adjusting);
            for (let i  = 0; i < tbs._data.length; i++) {
                tbs._data[i] = tbs._data[i].concat([0])
            }
            tbs = math.matrix(reducedREF(tbs))
            // EXTRACT BASIS FROM reduced ROW ECHELON FORM
            let ldothere = Array(tbs._size[0]).fill(false);
            let rer = 0;
            for (let rec = 0; rec < tbs._size[0]; rec++) {
                if (isEqualF(tbs.get([rer, rec]), 1)) {
                    ldothere[rec] = true;
                    rer++;
                }
                if (rer >= tbs.length) {
                    break;
                }
            }
            for (let rec = 0; rec < ldothere.length; rec++) {
                if (!ldothere[rec]) {
                    let toconc = math.transpose(math.column(tbs, rec))
                    toconc.set([0, rec], -1);
                    toconc = beautifyv(toconc);
                    tap.egvec = math.concat(tap.egvec, toconc)
                }
            }
            aoEP.push(tap);
        }
    }
    return aoEP;
}

function getSVD(m) {
    let theSVD = new SVDM();
    let leftsingular = math.multiply(m, math.transpose(m))
    let rightsingular = math.multiply(math.transpose(m), m)
    let lsEP = getEigenPair(leftsingular);
    let rsEP = getEigenPair(rightsingular);
    // U SECTION
    for (let i = 0; i < lsEP.length; i++) {
        let tap = normalizevect(math.row(lsEP[i].egvec, 0))
        if (theSVD.u._size[1] === 0) {
            theSVD.u = math.concat(theSVD.u, tap)
        } else {
            theSVD.u = math.concat(theSVD.u, tap, 0)
        }
    }
    theSVD.u = math.transpose(theSVD.u)
    // VT section
    for (let i = 0; i < rsEP.length; i++) {
        let tap = normalizevect(math.row(rsEP[i].egvec, 0))
        if (theSVD.vt._size[1] === 0) {
            theSVD.vt = math.concat(theSVD.vt, tap)
        } else {
            theSVD.vt = math.concat(theSVD.vt, tap, 0)
        }
    }
    // SIGMA SECTION
    theSVD.sigma = math.zeros(math.size(m));
    if (m._size[0] > m._size[1]) {
        // lsEP lebih gede
        for (let i = 0; i < theSVD.sigma._size[0]; i++) {
            if (!isEqualF(lsEP[i].egval, 0)) {
                theSVD.sigma.set([i, i], Math.sqrt(lsEP[i].egval));
            }
        }
    } else {
        for (let i = 0; i < theSVD.sigma._size[0]; i++) {
            if (!isEqualF(rsEP[i].egval, 0)) {
                theSVD.sigma.set([i, i], Math.sqrt(rsEP[i].egval));
            }
        }
    }
    return theSVD;
}

// console.log(mt)
// let mta = reducedREF(mt)
// console.log(mta)
// console.log(math.subtract(math.matrix([1,2,3]), math.matrix([2,3,4])))
// m = math.multiply(math.transpose(m), m)
// console.log(m)
// m.subset(math.index(1, [0,1,2]), [1,2,99])
// console.log(m)

