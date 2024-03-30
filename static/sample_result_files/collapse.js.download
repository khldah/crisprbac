$(document).ready(function() {
  $('.expandable').nextAll('tr').each(function() { // hide all at start
    if (!($(this).is('.expandable')))
      $(this).hide();
  });
  $('.expandable input[type=button]').click(function() { // toggle single by click
    var trElem = $(this).closest("tr");
    trElem.nextAll('tr').each(function() {
      if ($(this).is('.expandable')) {
        return false;
      }
      $(this).toggle();
    });
  });
  $('#expand_all').click(function() { // show all
    $('.expandable').nextAll('tr').each(function() {
      if (!($(this).is('.expandable')))
        $(this).show();
    });
  });
  $('#collaps_all').click(function() { // hide all
    $('.expandable').nextAll('tr').each(function() {
      if (!($(this).is('.expandable')))
        $(this).hide();
    });
  })
});
