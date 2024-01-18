window.onload = function(){
    let addReagentFieldBtn = document.getElementById('add-reagent-field');
    addReagentFieldBtn.addEventListener('click', update_table);

    let addProductFieldBtn = document.getElementById('add-product-field');
    addProductFieldBtn.addEventListener('click', update_table);

    const toggleFormBtn = document.getElementById('toggle-btn-form');
    toggleFormBtn.addEventListener('click', toggleFormHandler)
}


function update_table(e){
    let type_field;

    if(e.target.id === 'add-reagent-field'){
        type_field = 'reagents_info';
    }
    else{
        type_field = 'products_info';
    }

    e.preventDefault();
    let allTypeFieldWrapper = document.getElementById(type_field);
    let allTypeField = allTypeFieldWrapper.getElementsByTagName('input');

    if(allTypeField.length > 200) {
        alert('You can only add 200 items');
        return;
    }

    let ReagentInputIds = []

    for(let i = 0; i < allTypeField.length; i++) {
        ReagentInputIds.push(parseInt(allTypeField[i].name.split('-')[1]));

    }

    let newFieldName_name = `${type_field}-${Math.max(...ReagentInputIds) + 1}-reagent_name`;
    let newFieldName_role = `${type_field}-${Math.max(...ReagentInputIds) + 1}-reagent_role`;
    let newFieldName_mmass = `${type_field}-${Math.max(...ReagentInputIds) + 1}-molar_mass`;
    let newFieldName_mass = `${type_field}-${Math.max(...ReagentInputIds) + 1}-mass`;
    let newFieldName_cc = `${type_field}-${Math.max(...ReagentInputIds) + 1}-cc50`;
    let newFieldName_smiles = `${type_field}-${Math.max(...ReagentInputIds) + 1}-smiles`;

    allTypeFieldWrapper.insertAdjacentHTML('beforeend',`
    <table class='din_table'>
    <tr>
        <td><input class="form-control" id="${newFieldName_name}" name="${newFieldName_name}" required type="text" value=""></td>
        <td><input class="form-control" id="${newFieldName_role}" name="${newFieldName_role}" required type="text" value=""></td>
        <td><input class="form-control" id="${newFieldName_mmass}" name="${newFieldName_mmass}" required step="any" type="number" value=""></td>
        <td><input class="form-control" id="${newFieldName_mass}" name="${newFieldName_mass}" required step="any" type="number" value=""></td>
        <td><input class="form-control" id="${newFieldName_cc}" name="${newFieldName_cc}" step="any" type="number" value=""></td>
        <td><input class="form-control" id="${newFieldName_smiles}" name="${newFieldName_smiles}" type="text" value=""></td>
        <td><button class="btn btn-danger btn-sm">delete</button></td>
    </tr>
    </table>
    `);

    const deleteBtns = document.getElementsByClassName('btn btn-danger btn-sm');
    deleteBtns[deleteBtns.length - 1].addEventListener('click', (e) => {
        e.preventDefault();
        console.log(e.target.parentElement.parentElement.parentElement)
        e.target.parentElement.parentElement.parentElement.remove()
    })
}

function toggleFormHandler(e){
    const toggleFormBtn = document.getElementById('toggle-btn-form');
    const uploadFileForm = document.getElementById('upload-file-form');
    const createFileForm = document.getElementById('create-file-form');

    e.preventDefault()
    toggleFormBtn.classList.toggle('active')

    if (toggleFormBtn.classList.contains('active')) {
        toggleFormBtn.innerText = 'upload file'
        createFileForm.style.display = 'block'
        uploadFileForm.style.display = 'none'
    } else {
        toggleFormBtn.innerText = 'create file'
        createFileForm.style.display = 'none'
        uploadFileForm.style.display = 'block'
    }
}
