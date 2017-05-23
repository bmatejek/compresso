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
        var decompressed_data = COMPRESSO.Decompress(bytes);
    }
};