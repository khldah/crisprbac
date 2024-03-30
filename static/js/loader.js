document.onreadystatechange = function () {
    var state = document.readyState
    if (state == 'interactive') {
        document.getElementById('result').style.visibility="hidden";
    } else if (state == 'complete') {
        document.getElementById('interactive');
        document.getElementById('load').style.visibility="hidden";
        document.getElementById('result').style.visibility="visible";
    }
}