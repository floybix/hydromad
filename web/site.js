
/* Copyright (C) 2010 Felix Andrews <felix@nfrac.org> */

var currentAnchor = "INIT";
var currentItem = "intro";

// repeatedly check the anchor in the URL, to detect back/forward
function checkAnchor() {
    if (currentAnchor != document.location.hash) {  
        currentAnchor = document.location.hash;
	var newItem;
        if ((!currentAnchor) || (currentAnchor == "#")) {
	    newItem = "intro";
	} else {
            newItem = currentAnchor.substring(1);
        }
	loadItem(newItem);
    }
}  

function setAnchor(newItem) {
    var a;
    if (newItem == "intro") {
	a = "#";
    } else {
	a = "#" + newItem;
    }
    document.location.hash = a;
}

function loadItem(newItem) {
    // record in access log
//    pageTracker._trackPageview("/item/" + newItem);
    // animate change of page
    if (currentItem != newItem) {
	$(jq(currentItem)).slideUp();
	currentItem = newItem;
    }
    $(jq(newItem)).slideDown();
    // set menu item to 'active'
    $("#nav a.active").removeClass("active");
    var navEl = $("#nav " + jq("nav_" + newItem));
    navEl.addClass("active");
    if (newItem != "intro") {
	// expand the corresponding nav group
	openNavGroup(navEl.parents("li.navgroup"));
	// load man page immediately if there is no example image
	if ($(jq(newItem)).find("img").length == 0) {
	    $(jq(newItem)).find("a.helplink").click();
	}
    }
}

/* constructs an id selector, escaping '.' and ':' (from jquery.com) */
function jq(myid) { 
    return '#' + myid.replace(/(:|\.)/g,'\\$1');
}

function openNavGroup(el) {
    el.siblings("li.navgroup:visible").slideUp();
    el.slideDown();
    // set parent nav item to 'active'
    $("#nav li.navhead.active").removeClass("active");
    el.prev().addClass("active");
    el;
}

jQuery(function(){
	// suppress loading of images until they are needed:
	//$(".item img").attr("src", "");
	$(".item").hide();
	$(".groupname").hide();
	// collapse subnavigation initially
	$("#nav li.navgroup").hide();

	$("#nav li.navhead a").click(function() {
		openNavGroup($(this).parent().next());
		return false;
	    });

	$("#nav li.navgroup a").click(function(e) {
		e.preventDefault();
		var newItem = $(this).attr("id").substring(4);
		// set the URL anchor, will trigger loadItem()
		setAnchor(newItem);
	    });

	$("a.helplink").click(function(e) {
		e.preventDefault();
		helplink = $(this);
		href = helplink.attr("href");
		// record in access log
		//	pageTracker._trackPageview(href);
		helplink.slideUp();
		loader = $('<div class="loading">Loading...</div>');
		helplink.after(loader);
		loader.hide().slideDown();
		man = $('<div class="manpage"></div>');
		helplink.after(man);
		man.load(href, function() {
			loader.hide();
			$(this).find("h2,table:first,div:last").remove();
			$(this).hide().slideDown();
		    });
	    });

	checkAnchor();
	setInterval("checkAnchor()", 300);
});