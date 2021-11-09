let img_input = document.getElementById('input_img');
let file_input = document.getElementById('file_input');
const gpu = new GPU();
file_input.addEventListener('change', (e) => {
    img_input.src = URL.createObjectURL(e.target.files[0])
},false);

function multiplyMatrix(ma, mb) {
    let brow = mb.length;
    let resrow = ma.length;
    let rescol = mb[0].length;
    let res = gpu.createKernel(
        function(a, b, blen) {
            let sum = 0; 
            for (let i = 0; i < blen; i++) {
                sum += a[this.thread.y][i] * b[i][this.thread.x];
            }
            return sum;
        }
    ).setOutput([rescol, resrow]);
    let inresult = res(ma, mb, brow).map(i => Array.from(i));
    return inresult;
}
function arrToMat(arr) {
    let zeromat = cv.matFromArray(arr.length, arr[0].length, cv.CV_8UC1, [].concat(...arr))
    return zeromat;
}

function matToArr(mat) {
    let arr = [];
    let rows = mat.rows;
    let cols = mat.cols;
    for (let i = 0; i < rows; i++) {
        let rowNow = [];
        for (let j = 0; j < cols; j++) {
            rowNow.push(Array.from(mat.ucharPtr(i, j))[0]);
        }
        arr.push(rowNow);
    }
    return arr;
}

function compressChannel(channelArr, percent) {

    let {u, v, q} = svdoptimized(channelArr, 'f');
    let modifiedQ = [];
    let cols = channelArr[0].length;
    let level = Math.round(cols * percent / 100);
    for (let i = 0; i < level; i++) {
        let rowNow = [];
        for (let j = 0; j < level; j++) {
            if (i == j) {
                if (q[i] == NaN) {
                    rowNow.push(0);
                } else {
                    rowNow.push(q[i]);
                }
            } else {
                rowNow.push(0);
            }
        }
        modifiedQ.push(rowNow);
    }
    let vTranspose = math.transpose(math.matrix(v))._data;
    let modifiedVT = vTranspose.slice(0, level);
    let modifiedU = u.map(i => i.slice(0, level));
    let leftMatrix = multiplyMatrix(modifiedQ, modifiedVT);
    let finalArr = multiplyMatrix(modifiedU, leftMatrix).map(i => i.map(j => {
        if (j == NaN) {
            return 0;
        } else if (j > 255) {
            return 255;
        } else if (j < 0) {
            return -1*j;
        } else {
            return j;
        }
    }
    ));
    // let finalArr = math.multiply(math.matrix(u), math.multiply(math.matrix(modifiedQ), math.matrix(vTranspose)))._data;
    return finalArr;
}
// PROSES KOMPRESI GAMBAR
img_input.onload = function() {
    let src = cv.imread(img_input);
    let dst = new cv.Mat();
    let rgbaSrc = new cv.MatVector();
    let rgbaDst = new cv.MatVector();
    let zeroArr = [];
    for (let i = 0; i < src.rows; i++) {
        let rowNow = [];
        for (let j = 0; j < src.cols; j++) {
            rowNow.push(0);
        }
        zeroArr.push(rowNow);
    }
    let zeroMat = arrToMat(zeroArr);
    cv.split(src, rgbaSrc);
    for (let i = 0; i < 3; i++) {
        rgbaDst.push_back(arrToMat(compressChannel(matToArr(rgbaSrc.get(i)), 30)));
        console.log(matToArr(rgbaDst.get(i)));
        // rgbaDst.push_back(rgbaSrc.get(i));
    }
    // console.log(rgbaSrc.get(0));
    // console.log(rgbaSrc.get(1));
    // console.log(rgbaSrc.get(2));
    // rgbaDst.push_back(
    //     arrToMat(compressChannel(matToArr(rgbaSrc.get(0)), 1))
    // );
    // rgbaDst.push_back(zeroMat);
    // rgbaDst.push_back(zeroMat);
    // rgbaDst.push_back(rgbaSrc.get(3));
    // console.log(matToArr(rgbaDst.get(0)));
    // let RMat = rgbaSrc.get(0);
    // let GMat = rgbaSrc.get(1);
    // let BMat = rgbaSrc.get(2);
    // let RArr = matToArr(RMat);
    // let GArr = matToArr(GMat);
    // let BArr = matToArr(BMat);
    // // let AArr = matToArr(AMat);
    // let complevel = 150
    // let RArrComp = compressChannel(RArr, complevel);
    // let GArrComp = compressChannel(GArr, complevel);
    // let BArrComp = compressChannel(BArr, complevel);
    // // let AArrComp = compressChannel(AArr, complevel);
    // let RMatSec = arrToMat(RArrComp);
    // let GMatSec = arrToMat(GArrComp);
    // let BMatSec = arrToMat(BArrComp);
    // // let AMatSec = arrToMat(AArrComp);
    // rgbaDst.push_back(RMatSec);
    // rgbaDst.push_back(GMatSec);
    // rgbaDst.push_back(BMatSec);
    // rgbaDst.push_back(AMatSec);
    cv.merge(rgbaDst, dst);
    cv.imshow('output', dst);
    src.delete();
}