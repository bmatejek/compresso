var SEGMENTATION = SEGMENTATION || {};

SEGMENTATION.File = function(file) {
    this._file = file;

    /* create a reader variable */
    this._reader = new FileReader();
    this._reader.onerror = this.onError.bind(this);
    this._reader.onloadend = this.onLoad.bind(this);

    /* read the entire file */
    this._reader.readAsArrayBuffer(file);
};



SEGMENTATION.File.prototype.onError = function(e) {
    throw new Error(e.target.error.code);
};



SEGMENTATION.File.prototype.onLoad = function(e) {
    if (e.target.readyState == FileReader.DONE) {
        var buffer = e.target.result;

        var bytes = new Uint32Array(buffer);

        // call COMPRESSO algorithm
        var data = COMPRESSO.Decompress(bytes);

        var canvas = document.createElement("canvas");
        canvas.width = data.xres;
        canvas.height = data.yres;
        var ctx = canvas.getContext("2d");

        var iz = 4;
        for (var ix = 0; ix < data.xres; ++ix) {
            for (var iy = 0; iy < data.yres; ++iy) {
                component = data.decompressed_data[IndicesToIndex(ix, iy, iz)];

                red = (Math.trunc(107 * component) % 700) % 255;
                green = (Math.trunc(509 * component) % 900) % 255;
                blue = (Math.trunc(200 * component) % 777) % 255;

                ctx.fillStyle = 'rgb(' + red + ', ' + green + ', ' + blue + ')';
                ctx.fillRect(ix, iy, ix + 1, iy + 1);
            }
        }

        var img = document.createElement("img");
        img.src = canvas.toDataURL("image/png");
        document.body.appendChild(img);
    }
};