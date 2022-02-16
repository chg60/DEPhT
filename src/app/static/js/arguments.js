/**Show optional arguments.*/
function show_optional(opt_id, show_id, hide_id) {
  opt = document.getElementById(opt_id);
  show = document.getElementById(show_id)
  hide = document.getElementById(hide_id)

  if(opt.style.display === "none" || opt.style.display == "") {
    opt.style.display = "block";
    show.style.display = "none";
    hide.style.display = "block";
  } else {
    opt.style.display = "none";
    show.style.display = "block";
    hide.style.display = "none";
  }
}


/**
* Error message for required arguments.
*
* @param {Array} id_list
*/
function validity_check(id) {
  element = document.getElementById(id);
  validity_state = element.validity;

  window.alert(id);

  if(validity_state.valueMissing) {
    element.setCustomValidity("This field is required.");
    element.style.border = "solid red"
  } else if(validity_state.rangeUnderflow || validity_state.rangeOverflow) {
    element.setCustomValidity("Out of range.");
  } else {
    element.setCustomValidity("");
  }

  element.reportValidity();
}


/**
* Get list of arguments.
*
*/
function get_args_list() {
  args_list = document.querySelectorAll("select, input");
  const fs = require("fs");

  for(let i = 0; i < args_list.length; i++) {
    // .type, .id, .value
    fs.writeFile("args.txt", args.id+","+args.value, (err) => {
      if (err) throw err;
    });
  }

}

function check_input_validity() {
  id_list = document.querySelectorAll("select, input");

  id_list.forEach(id => {
    validity_check(id);
  });
}


/**
* Get a list of number of cores.
*/
function get_core_list(cpu_cores) {
  select = document.getElementById("cores");

  for(let i = 1; i <= cpu_cores; i++) {
    var option = document.createElement("option");
    option.text = option.value = i;
    select.add(option);
  }
}


function loading_page() {
  window.location.href = "/results";
}
