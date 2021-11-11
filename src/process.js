const userDataForm = document.querySelector('#compressionsettings');
const comprateinp = document.querySelector('#compRateinp');
const imageInp = document.querySelector('#imageinp');
const before_img = document.querySelector('#before_img');
const isTruePixel = document.querySelector('#isTruePixel');
const truepx_img = document.querySelector('#truePxImg');
const downloadbutton = document.querySelector('.dlbtn');
const perfstats = document.querySelector('.perfstats');
let filename;
let filetype;
let startTime;
let endTime;
let aftercanvas = document.querySelector('#after_img');
aftercanvas.width = 0;
aftercanvas.height = 0;
let origwidth = 0;
let origheight = 0;
const gpu = new GPU();
let beforeSrc = '';
prevSvdMatrices = [];

isTruePixel.addEventListener('change', function(e) {
    prevSvdMatrices = [];
})

function dlCanvas() {
    let linkdl = document.createElement('a');
    let origfilename = filename.substring(0, filename.lastIndexOf("."));
    let origext = filename.substring(filename.lastIndexOf(".") + 1, filename.length);
    linkdl.download = `${origfilename}_compressed.${origext}`;
    linkdl.href = aftercanvas.toDataURL(filetype);
    linkdl.click();
}

function dlCanvasResized() {
    let resizedcv = document.createElement('canvas');
    let resizedctx = resizedcv.getContext('2d');
    resizedcv.width = origwidth;
    resizedcv.height = origheight;
    resizedctx.drawImage(aftercanvas, 0, 0, origwidth, origheight);
    let linkdl = document.createElement('a');
    let origfilename = filename.substring(0, filename.lastIndexOf("."));
    let origext = filename.substring(filename.lastIndexOf(".") + 1, filename.length);
    linkdl.download = `${origfilename}_compressed.${origext}`;
    linkdl.href = resizedcv.toDataURL(filetype);
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

function displayPerfStats() {
    perfstats.innerHTML = `
    <hr/>
    <p class="card-text" style="margin-bottom: 5px !important;"><small class="text-muted" id="performancetimestats">Time: ${Math.round((Math.round(endTime - startTime) / 1000) * 100)/100}s</small>
    </p>
    <p class="card-text"><small class="text-muted" id="compratestats">Compression Rate: ${101 - comprateinp.value}%</small></p>
    `
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
    let zeroatfirst = 0;
    for (let i = 0; i < q.length; i++) {
        if (q[i] > 0) {
            break;
        } else {
            zeroatfirst++;
        }
    }
    let level = rows < cols ? Math.round((rows - zeroatfirst) * percent / 100) : Math.round((cols - zeroatfirst) * percent / 100);
    level = level + zeroatfirst;
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
    filename = e.target.files[0].name;
    filetype = e.target.files[0].type;
    beforeSrc = URL.createObjectURL(e.target.files[0]);
    prevSvdMatrices = [];
},false);

userDataForm.addEventListener('submit', function(e) {
    e.preventDefault();
    if (isTruePixel.checked) {
        truepx_img.cprate = comprateinp.value;
        truepx_img.src = beforeSrc;
        before_img.src = beforeSrc;
    } else {
        before_img.src = beforeSrc;
        before_img.cprate = comprateinp.value;
    }
}); 

before_img.onload = function() {
    if (!isTruePixel.checked) {
        origwidth = this.naturalWidth;
        origheight = this.naturalHeight;
        startTime = performance.now();
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
        endTime = performance.now();
        src.delete();
        dst.delete();
        displayPerfStats();
        genDlButton();
    }
}

truepx_img.onload = function() {
    if (isTruePixel.checked) {
        startTime = performance.now();
        let src = cv.imread(truepx_img);
        let dst = new cv.Mat();
        let rgbaSrc = new cv.MatVector();
        let rgbaDst = new cv.MatVector();
        cv.split(src, rgbaSrc);
        for (let i = 0; i < 3; i++) {
            let out = compressChannel(matToArr(rgbaSrc.get(i)), truepx_img.cprate, i);
            rgbaDst.push_back(arrToMat(out));
        }
        rgbaDst.push_back(rgbaSrc.get(3))
        cv.merge(rgbaDst, dst);
        cv.imshow('after_img', dst);
        endTime = performance.now();
        src.delete();
        dst.delete();
        displayPerfStats();
        genDlButtonTruePixel();
    }
}