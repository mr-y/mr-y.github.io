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
    //var n_branches = 0;
    this.scores = {};
    this.copy = function () {
	treeCopy = new tree();
	treeCopy.root = copyComplexObject(this.root);
	treeCopy.scores = copyComplexObject(this.scores);
    }
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
		//++n_branches;
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
    this.write = function() {
	var output = "";
	function writeNode(node) {
	    if (node.children.length > 0) { output += '('; }
	    for (var i=0; i < node.children.length; ++i) {
		if (i !== 0) { output += ','; }
		writeNode(node.children[i]);
	    }
	    if (node.children.length > 0) { output += ')'; }
	    output += node.name;
    	    output += ":" + node.branch_length;
	}
	writeNode(this.root);
	output += ';';
	return output;
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
    this.calc_parsimony_scores = function (alignment, alphabet) {
	//if (this.scores === undefined) { this.scores = {}; }
	var translation = alphabet;
	this.add_parsimony = function (node) {
	    if (this.scores === undefined) { alert("this.scores === undefined"); }
	    if (node.children.length > 0) {
    		for (var i=0; i < node.children.length; ++i) {
		    this.add_parsimony(node.children[i]);
		}
	    }
	    if (node.name && alignment.OTUs[node.name]) {
		node.seq = alignment.OTUs[node.name];
		for (var i=0; i < node.children.length; ++i) {
		    var done_part = {};
                    for (var part in node.seq) {
			if (node.seq.hasOwnProperty(part) && node.children[i]['seq'].hasOwnProperty(part)) {
			    if (this.scores[part] === undefined) { this.scores[part] =[]; }
			    for (var pos=0; pos < node.seq[part].length; ++pos) {
				var comp = translation.pars_compare(node.seq[part][pos],node.children[i]['seq'][part][pos]);
			       	if (comp[1]) {
				    if (!this.scores[part][pos]) { this.scores[part][pos] = 0; }
				    this.scores[part][pos] += comp[1];
			       	}
			    }
			    done_pars[part] = true;
			}
		    }
		    for (var part in node.children[i]['seq']) {
			if (node.children[i]['seq'].hasOwnProperty(part) && !done_part[part]) {
			    node.seq[part] = [];
			    for (var pos=0; pos < node.children[i]['seq'][part].length; ++pos) {
				node.seq[part][pos] = node.children[i]['seq'][part][pos];
			    }
			}
		    }
                }
	    }
	    else {
		if (!node.seq) { node.seq = {}; }
		for (var i=0; i < node.children.length; ++i) {
		    for (var part in node.children[i]['seq']) {
			if (node.children[i]['seq'].hasOwnProperty(part)) {
			    if (this.scores[part] === undefined) { this.scores[part] = []; }
			    if (node.seq[part]) {
				for (var pos=0; pos < node.children[i]['seq'][part].length; ++pos) {
				    var comp = translation.pars_compare(node.seq[part][pos],node.children[i]['seq'][part][pos]);
				    //alert (node.seq[part][pos] + " " + node.children[i]['seq'][part][pos] + " " + comp);
				    node.seq[part][pos] = comp[0];
				    if (!this.scores[part][pos]) { this.scores[part][pos] = 0; }
				    this.scores[part][pos]+= comp[1];
				}
			    }
			    else {
				node.seq[part] = [];
				for (var pos=0; pos < node.children[i]['seq'][part].length; ++pos) {
				    node.seq[part][pos] = node.children[i]['seq'][part][pos];
				}
			    }
			}
		    }
		}
	    }
	}
	this.add_parsimony (this.root);
    }
    this.tot_score = function () {
	var sum = 0;
	for (var part in this.scores) {
    	    if (this.scores.hasOwnProperty(part)) {
		for (var i=0; i < this.scores[part].length; ++i) {
		    sum += this.scores[part][i];
		}
	    }
	}
	return sum;
    }
    this.nni = function (branch) {
	var newTrees = [];
	var newBranch=new node();
	//console.log("Branch: " + branch);
	function copyWithPointer(A, pointer, mother=null) {
	    if (Array.isArray(A)) {
                //console.log("Copy array");
                var B = [];
                for (var i=0; i< A.length; ++i) {
                    if (A[i] instanceof Array || A[i] instanceof Object) {
                        B[i] = copyWithPointer(A[i], pointer, mother);
                    }
                    else { B[i] = A[i]; }
                }
                return B;
            }
	    else if (A instanceof Object) {
		var B = {};
		if (A.mother) { B.mother = mother; mother = B; }
		for (part in A) {
		    if (A.hasOwnProperty(part) && part !== 'mother') {
			if (A[part] instanceof Array || A[part] instanceof Object) {
			    //console.log(part + ' ' + typeof A[part]);
			    B[part] = copyWithPointer(A[part], pointer, mother);
			}
			else { B[part] = A[part]; }
		    }
		}
		if (A === pointer) { newBranch = B; }//console.log("Found branch"); }
		//console.log("New pointer children: " + newBranch.children)
		return B;
	    }
	}
	if (!branch || branch === this.root || branch.children.length < 1) { console.log("Root"); return newTrees; }
	if (branch.mother !== undefined && branch.mother !== null && branch.mother !== 0 && branch.children.length > 1) {
	    if (branch.mother !== this.root) {
		newTrees[0] = new tree();
		newTrees[0].scores = {};
	       	newTrees[0].root = copyWithPointer(this.root,branch);
		//console.log('Test: ' + Array.isArray(this.root.mother) );
		var tempNode = new node();
		/*for (part in newBranch) {
		    if (newBranch.hasOwnProperty(part)) console.log(part + ' ' + newBranch[part]);
		}*/
		tempNode = newBranch.children[1];
		//console.log("YY: " + tempNode.children);
		newBranch.children[1] = newBranch.children[0];
		if (newBranch.mother.children[0] === newBranch) { newBranch.children[0] = newBranch.mother.children[1]; newBranch.mother.children[1] = tempNode; }
		else { newBranch.children[0] = newBranch.mother.children[0]; newBranch.mother.children[0] = tempNode; }
		newTrees[1] = new tree();
		newTrees[1].scores = {};
	       	newTrees[1].root = copyWithPointer(this.root,branch);
		tempNode = newBranch.children[0];
		newBranch.children[0] = newBranch.children[1];
		if (newBranch.mother.children[0] === newBranch) { newBranch.children[1] = newBranch.mother.children[1]; newBranch.mother.children[1] = tempNode; }
		else { newBranch.children[1] = newBranch.mother.children[0]; newBranch.mother.children[0] = tempNode; }
		
	    }
	    else {
		//console.log("First child: " + this.root.children[0]);
		var rootChild;
		if (this.root.children[0] === branch) {
		    rootChild=1;
		    console.log("Left");
		    if (this.root.children[1].children.length < 2) { return newTrees; }
		}
		else if (this.root.children[1] === branch) {
		    rootChild=0;
		    console.log("Right");
		    if (this.root.children[0].children.length < 2) { return newTrees; }
		}
		console.log("Gren " + rootChild);
		newTrees[0] = new tree();
	    	newTrees[0].scores = {};
		newTrees[0].root = copyWithPointer(this.root,branch);
		tempNode = new node();
		tempNode = newBranch.children[1];
		newBranch.children[1] = newTrees[0].root.children[rootChild].children[0];
		root.children[rootChild].children[0] = tempNode;
		newTrees[1] = new tree();
		newTrees[1].scores = {};
		newTrees[1].root = copyWithPointer(this.root,branch);
		tempNode = newBranch.children[1];
		newBranch.children[1] = newBranch.children[0].children[1];
		root.children[1].children[1] = tempNode;
	    }
	}
	return newTrees;
    }

    this.getNodesAsArray = function () {
	var nodes = [];
	function getNodes (node) {
	    for (var i = 0; i < node.children.length; ++i) {
		getNodes(node.children[i]);
	    }
	    nodes.push(node);
	}
	getNodes(this.root);
	return nodes;
    }

    /*this.nni_all = function () {
	var nniTrees = []
	function  transverse (node) {
	    for (var i = 0; i < node.children.length; ++i) {
		this.transverse(node.children[i]);
	    }
	    nniTrees.push(this.nni(node));
	}
	this.transverse(this.root);
	return nniTrees;
    }*/

    this.add_svg_annotation = function (width,height) {
	var nTaxa = this.nTips();
	var perTaxa = (height/*-(height/(2*nTaxa))*/)/nTaxa;
	var xScale = (width-width*0.1)/this.height();
	var line_width = 3;
	function sub_add_svg (node, x_start, toTheLeft) {
	    var y_start;
	    if (node.children.length > 0) {
		var xy_end = [];
		for (var i=0; i<node.children.length;++i){
		    xy_end[i] = sub_add_svg(node.children[i], x_start + node.branch_length*xScale,toTheLeft);
		}
		if (xy_end) { y_start = (xy_end[0][1]+xy_end[xy_end.length-1][1])/2; }
		x_start = x_start+node.branch_length*xScale;
		for (var i=0; i<node.children.length;++i) {
		    if (!node.children[i]['svg']) node.children[i]['svg'] = {};
		    node.children[i]['svg']['polyline'] = [];
		    node.children[i]['svg']['polyline'][0] = {};
		    node.children[i]['svg']['polyline'][0]['atributes'] = "points=\"" + x_start + "," + y_start + " " + x_start + "," + xy_end[i][1] + " " + xy_end[i][0] + "," + xy_end[i][1] + "\"";
		    node.children[i]['svg']['polyline'][0]['atributes'] += " style=\"fill:none;stroke:black;stroke-width:" + line_width + "\"";
		    node.children[i]['svg']['polyline'][0]['atributes'] += " id=\"node_" + toTheLeft.nTips + "\"";
		}
	    }
	    else {
		x_start = x_start+node.branch_length*xScale; y_start = (toTheLeft.nTips*perTaxa)+(0.5*perTaxa);
		++toTheLeft.nTips;
	    }
	    if (node.name) {
		if (!node['svg']) node['svg'] = {};
		node['svg']['text'] = [];
		node['svg']['text'][0] = {};
		node['svg']['text'][0]['atributes'] = "x=\"" + x_start + "\" y=\"" + (y_start+5) + "\"";
		node['svg']['text'][0]['value'] = node.name;
	    }
    	    return [x_start,y_start];
	}
	sub_add_svg(this.root,0,{nTips:0});
    }
    this.drawSVG = function ( width,height) {
	var SVGdrawing = "<svg height=\"" + height + "\" width=\"" + width + "\">\n";
	function subDraw (node) {
	    if (node.hasOwnProperty('svg')) {
		for (element in node['svg']) {
		    if (node['svg'].hasOwnProperty(element)) {
			for (var i=0; i < node['svg'][element].length; ++i) {
			    SVGdrawing += "<" + element;
			    if (node['svg'][element][i]['atributes']) SVGdrawing += " " + node['svg'][element][i]['atributes'];
			    SVGdrawing += ">";
			    if (node['svg'][element][i]['value']) SVGdrawing += node['svg'][element][i]['value'];
			    SVGdrawing += "</" + element + ">";
			}
		    }
		}
	    }
	    for (var i=0; i < node.children.length; ++i) {
		subDraw(node.children[i]);
	    }
	}
	subDraw (this.root);
	SVGdrawing += "</svg>\n";
	return SVGdrawing;
    }
}

function copyComplexObject(A,mother=null) {
    var B = {};
    for (part in A) {
	if (A.hasOwnProperty(part)) {
	    if (part === 'mother') { B[part] = mother; }
	    else if (A[part] instanceof Array) {
		console.log("Object " + part );
		B[part] = copyComplexArray(A[part], B);
	    }
	    else if (A[part] instanceof Object) {
		console.log("Array " + part);
		B[part] = copyComplexObject(A[part], B);
	    }
	    else { B[part] = A[part]; console.log(B[part] + " " +part);}
	}
    }
    return B;
}

function copyComplexArray(A,mother=null) {
    var B = [];
    for (var i=0; i < A.length; ++i) {
	if (A[i] instanceof Array) {
	    console.log("Array " + part);
	    B[i] = copyComplexArray(A[i]);
	}
	else if (A[i] instanceof Object) {
	    console.log("Object " + part);
	    B[i] = copyComplexObject(A[i], mother);
	}
	else { B[i] = A[i]; console.log(B[i] + " " +part);}
    }
    return B;
}

