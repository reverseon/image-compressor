const userDataForm = document.querySelector('#compressionsettings');
const comprateinp = document.querySelector('#compRateinp');
const imageInp = document.querySelector('#imageinp');
const before_img = document.querySelector('#before_img');
const isTruePixel = document.querySelector('#isTruePixel');
const truepx_img = document.querySelector('#truePxImg');
const downloadbutton = document.querySelector('.dlbtn');
let aftercanvas = document.querySelector('#after_img');
aftercanvas.width = 0;
aftercanvas.height = 0;
let origwidth = 0;
let origheight = 0;
const gpu = new GPU();
let beforeSrc = '';
prevSvdMatrices = [];

function dlCanvas() {
    let linkdl = document.createElement('a');
    linkdl.download = 'compressed.png';
    linkdl.href = aftercanvas.toDataURL('image/png');
    linkdl.click();
}

function dlCanvasResized() {
    let resizedcv = document.createElement('canvas');
    let resizedctx = resizedcv.getContext('2d');
    resizedcv.width = origwidth;
    resizedcv.height = origheight;
    resizedctx.drawImage(aftercanvas, 0, 0, origwidth, origheight);
    let linkdl = document.createElement('a');
    linkdl.download = 'compressed.png';
    linkdl.href = resizedcv.toDataURL('image/png');
    linkdl.click();
}

function genDlButtonTruePixel() {
    downloadbutton.innerHTML = `
    <hr/>
    <button type="button" id="downloadcanvas" class="btn btn-primary float-end">Download</button>
    `
}

function genDlButton() {
    downloadbutton.innerHTML = `
    <hr/>
    <button type="button" id="downloadcanvas" class="btn btn-success float-start">Download</button>
    <button type="button" id="downloadcanvasresized" class="btn btn-primary float-end">Download Resized</button>
    `;
}

downloadbutton.addEventListener('click', function(e) {
    if (e.target.id == 'downloadcanvas') {
        dlCanvas();
    } else if (e.target.id == 'downloadcanvasresized') {
        dlCanvasResized();
    }
});

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

function compressChannel(channelArr, percent, i) {
    if (percent < 1 || percent > 100) {
        throw new Error('Percent must be between 1 and 100');
    }
    if(prevSvdMatrices.length < i+1){
        let svd = svdgolub(channelArr);
        prevSvdMatrices.push(svd);
        return multiplySvdToMat(svd, percent);
    } else{
        return multiplySvdToMat(prevSvdMatrices[i], percent);
    }
}

function multiplySvdToMat(svd, percent){
    let modifiedQ = [];
    let {u, q, v} = svd;
    let cols = svd.u.length;
    let rows = svd.v.length;
    let level = rows < cols ? Math.round(rows * percent / 100) : Math.round(cols * percent / 100);

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
    let vTranspose = math.transpose(v);
    let modifiedVT = vTranspose.slice(0, level);
    let modifiedU = u.map(i => i.slice(0, level));
    let leftMatrix = multiplyMatrix(modifiedQ, modifiedVT);
    let finalArr = multiplyMatrix(modifiedU, leftMatrix).map(i => i.map(j => {
        if (j == NaN || j == Infinity) {
            return 255;
        } else if (j > 255) {
            return 255;
        } else if (j < 0) {
            return -1*j;
        } else {
            return j;
        }
    }
    ));
    return finalArr;
}

imageInp.addEventListener('change', (e) => {
    beforeSrc = URL.createObjectURL(e.target.files[0]);
    prevSvdMatrices = [];
},false);

userDataForm.addEventListener('submit', function(e) {
    e.preventDefault();
    before_img.src = beforeSrc;
    before_img.cprate = comprateinp.value;
    if (isTruePixel.checked) {
        truepx_img.cprate = comprateinp.value;
        truepx_img.src = beforeSrc;
    }
}); 

before_img.onload = function() {
    origwidth = this.naturalWidth;
    origheight = this.naturalHeight;
    if (!isTruePixel.checked) {
        let src = cv.imread(before_img);
        let dst = new cv.Mat();
        let rgbaSrc = new cv.MatVector();
        let rgbaDst = new cv.MatVector();
        cv.split(src, rgbaSrc);
        for (let i = 0; i < 3; i++) {
            rgbaDst.push_back(arrToMat(compressChannel(matToArr(rgbaSrc.get(i)), before_img.cprate, i)));
        }
        rgbaDst.push_back(rgbaSrc.get(3))
        cv.merge(rgbaDst, dst);
        cv.imshow('after_img', dst);
        src.delete();
        dst.delete();
        genDlButton();
    }
}

truepx_img.onload = function() {
    if (isTruePixel.checked) {
        let src = cv.imread(truepx_img);
        let dst = new cv.Mat();
        let rgbaSrc = new cv.MatVector();
        let rgbaDst = new cv.MatVector();
        cv.split(src, rgbaSrc);
        for (let i = 0; i < 3; i++) {
            rgbaDst.push_back(arrToMat(compressChannel(matToArr(rgbaSrc.get(i)), truepx_img.cprate, i)));
        }
        rgbaDst.push_back(rgbaSrc.get(3))
        cv.merge(rgbaDst, dst);
        cv.imshow('after_img', dst);
        src.delete();
        dst.delete();
        genDlButtonTruePixel();
    }
}