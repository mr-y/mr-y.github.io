function node (label="", length=0, mom=null) {
    this.children=[];
    this.mother = mom;
    this.name = label;
    this.branch_length = length;
}

function tree () {
    this.root = new node();
    //this.node = new node();
    //this._root = node;
    var n_branches = 0;
    this.pars_newick = function (tree) { // tree should be a textstring witha newick formated tree
	//var present = this._root;
	var present = this.root;
	var read_mode = 's';
	var label="";
	var branch_length="";
	for (var i=0; i < tree.length; ++i) {
	    if (tree[i] === '[') {
		var n_square_right = 0;
		var start=true;
		while ((n_square || start) && i < tree.length) {
		    if (start) { start=false; }
		    if (tree[i] === '[') { ++n_square_right; }
		    else if (tree[i] === ']') { --n_square_right; }
		    ++i;
		}
	    }
	    else if (tree[i] === '(') {
		read_mode = 'l';
		present.children[present.children.length] = new node("",0,present);
		present = present.children[present.children.length-1];
		label="";
		branch_length="";
	    }
	    else if (tree[i] === ',') {
		++n_branches;
		read_mode = 'l';
	       	if (label.length > 0) {
		    present.name=label;
		}
		if (branch_length.length > 0) {
		    present.branch_length = parseFloat(branch_length);
		}
		label="";
		branch_length="";
		if (present.mother) { present = present.mother; }
		else { alert ("Failure to pars tree at char: " + i + " (" + tree[i] + ")" ) }
		present.children[present.children.length] = new node(0,0,present);
		present = present.children[present.children.length-1];
    	    }
    	    else if (tree[i] === ')') {
		read_mode = 'l';
		if (label.length > 0) {
                    present.name=label;
                }
                if (branch_length.length > 0) {
                    present.branch_length = parseFloat(branch_length);
                }
                label="";
                branch_length="";
		present = present.mother;
    	    }
    	    else if (tree[i] === ':') {
		read_mode = 'b';
    	    }
    	    else if (tree[i] === ';') {
		if (present === this.root) {
		//if (present === this._root) {
		    if (label.length > 0) {
			present.name=label;
		    }
		    if (branch_length.length > 0) {
			present.branch_length = parseFloat(branch_length);
		    }
		}
		else { alert ("Tree may not have been parsed correctly. Check tree format."); }
    	    }
    	    else if (read_mode === 'b') {
		branch_length += tree[i];
    	    }
    	    else if (read_mode === 'l') {
		if (label.length > 0 && (tree[i] === "'" || tree[i] === '"')) {
		    var marker = tree[i];
		    ++i;
		    while (tree[i] !== marker && i < tree.length) {
			label += tree[i];
		    }
		}
		else { label += tree[i]; }
    	    }
	}
    }
    this.nTips = function () {
	function subNtips(node) {
	    var n=0;
	    if (node.children.length < 1) { ++n; }
	    else {
		for (var i=0; i < node.children.length; ++i) {
		     n += subNtips(node.children[i])
		}
	    }
	    return n;
	}
	var n=0;
	return subNtips(this.root);
    }
    this.nNodes = function () {
	function subNnodes(node) {
	    var n=1;
	    for (var i=0; i < node.children.length; ++i) {
   		n += subNnodes(node.children[i])
	    }
	    return n;
	}
	var n=0;
	return subNnodes(this.root);
    }
    this.length = function () {
	function subLength (node) {
	    var length = 0;
	    for (var i=0; i < node.children.length; ++i) {
                     length += subLength(node.children[i])
	    }
	    return length+node.branch_length;
	}
	return subLength(this.root);
    }
    this.height = function () {
	function subHeight (node) {
	    var height  = 0;
	    for (var i=0; i < node.children.length; ++i) {
	       	var length = subHeight(node.children[i]);
		if (length > height) { height = length; }    
	    }
	    return height+node.branch_length;
	}
	return subHeight(this.root);
    }
    this.colless = function () {
	function subColless (node) {
	    var values = [];
	    var sum = [0, 0];
	    for (var i=0; i < node.children.length; ++i) {
		values[i] = subColless(node.children[i]);
		sum[0] += values[i][0];
		sum[1] += values[i][1];
	    }
	    if (node.children.length === 2) {
		sum[0] += Math.abs(values[0][1]-values[1][1]);
	    }
	    ++sum[1];
	    return sum;
	}
	return subColless(this.root)[0];
    }

    this.B1 = function () {
	function subB1 ( node ) {
	    var values = [];
	    var return_values = [0,0];
	    for (var i=0; i < node.children.length; ++i) {
		values = subB1(node.children[i]);
		if (values[1] > return_values[1]) { return_values[1] = values[1]; }
		return_values[0] += values[0];
	    }
	    if (node.children.length > 1) { return_values[0] += 1/return_values[1]; }
	    ++return_values[1];
	    return return_values;
	}
	var value = 0;
	for (var i=0; i < this.root.children.length; ++i) {
	    value += subB1(this.root.children[i])[0];
	}
	return value;
    }
    this.node_depths = function () {
	var nodedepths = [];
	function getNodeDepths (node, length) {
	    var depth = length+node.branch_length;
	    if (node.children.length > 1) { nodedepths.push(depth); }
	    for (var i=0; i < node.children.length; ++i) {
		getNodeDepths(node.children[i], depth);
	    }
	}
	getNodeDepths(this.root, 0);
	nodedepths.push(this.height());
	nodedepths.sort(function(a, b){return a-b});
	return nodedepths;
    }
    this.gamma = function () {
	var nodedepths = this.node_depths ();
	var sum = 0;
	var T = 0;
	for (var i=1; i < nodedepths.length; ++i) {
	    T += (i+1) * (nodedepths[i] - nodedepths[i-1]);
	    if (i < nodedepths.length-1) {
		for (var k=1; k<=i; ++k) {
		    sum += (k+1) * (nodedepths[k] - nodedepths[k-1]);
		}
	    }
	}
	var value = (1/(nodedepths.length-2)) * sum;
	value -= T/2;
	value /= T * Math.sqrt(1/(12*(nodedepths.length-2)));
	return value;
	//return nodedepths;
    }
    /*this.add_svg_annotation = function () {
	function sub_add_svg (node, n_left, scale_w, scale_x) {
	    var n_tot;
	    var daughters;
	    for (var i=0;;) {}
	}
    }*/
}
/*
function svg_annotations () {
    function polyline (line) {
	this.points = line;
	this.style;
	this.add_to_drawing = function (drawing) {
	    drawing += "<polyline points=\"" + this.points + "\"";
	    if (this.style) drawing += " style=\"" + this.style + "\"";
	    drawing += " />\n";
	}
    }
    function text (text, x, y) {
	this.text = text;
	this.x = x;
	this.y = y;
	this.fill;
	this.transform;
	this.add_to_drawing = function (drawing) {
	    drawing += "<text x=\"" + x + "\" y=\"";
	    if (this.fill) drawing += " fill=\"" + this.fill + "\"";
	    if (this.transform) rawing += " transform=\"" + this.transform + "\"";
	    drawing += ">" + this.text + "</text>";
	}
    }
    this.objects = [];
    this.add_branch = function (branch_length, start, scale=1, end) {
	this.objects.push(new polyline(start[0] + "," + start[1] + " " + start[0] + "," + (start[0]- branch_length * scale) + end[0] + "," + end[1]))
    }
    this.add_branch_label = function (label, branch_length, start, scale=1, x_offset, y_offset) {
	this.objects.push(new text (label, start[0]-x_offset, start[1] + branch_length * scale - y_offset ));
    }
    this.add_tip_label = function(label, start, x_offset, y_offset) {
	this.objects.push(new text (label, start[0]-x_offset, start[1] - y_offset ));
    }
}

function parsimony () {
    this.tree = new tree();
    this.alignment = new alignment();
    this.score_partitions = function () {
	

    }
}
*/
