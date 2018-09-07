firebase.initializeApp(FIREBASE_CONFIG);

firebase.firestore().settings({
    // Disable deprecated features
    timestampsInSnapshots: true
});


firebase.firestore().enablePersistence()
    .catch(function(err) {
        if (err.code == 'failed-precondition') {
            // Multiple tabs open, persistence can only be enabled
            // in one tab at a a time.
            // ...
        } else if (err.code == 'unimplemented') {
            // The current browser does not support all of the
            // features required to enable persistence
            // ...
        }
    });


const db = firebase.firestore();
const ref = db.collection("gtf_entries");


const geneHeader = document.querySelector('#gene');
const annotDiv = document.querySelector('#annot');

const inputTextField = document.querySelector('#currentGeneInput');
const saveButton = document.querySelector('#saveButton');

console.log(geneHeader, annotDiv, saveButton);


saveButton.addEventListener("click", () => {
    const geneName = inputTextField.value;
    console.log(`fetching annotation for ${geneName}`);

    let beg = performance.now();
    let query = ref.where('gene_name', '==', geneName).orderBy('start');
    query.get()
        .then((r) => {
            let end = performance.now();
            geneHeader.innerText = geneName + `(${end - beg}ms)`;
            let entries = r.docs.map((d) => d.data()).filter((d) => d.feature != "gene" & d.feature != "transcript");
            let html = entries.map(
                (e, i) => `<tr>` +
                    `<td>${i + 1}` +
                    `<td>${e.feature}` +
                    `<td>${e.start}</td>` +
                    `<td>${e.end}</td>` +
                    `<td>${e.transcript_id}</td>` +
                    `<td>${e.gene_name}</td>` +
                    `<td>${e.gene_id}</td>` +
                    `<td>${e.strand}</td>` +
                    `<td>${e.seqname}</td>` +
                    `</tr>`
            ).join(' ');

            let header = `<tr>` +
                `<th>#</th>` +
                `<th>feature</th>` +
                `<th>start</th>` +
                `<th>end</th>` +
                `<th>transcript id</th>` +
                `<th>gene name</th>` +
                `<th>gene id</th>` +
                `<th>strand</th>` +
                `<th>seqname</th>` +
                `</tr>`;

            annotDiv.innerHTML = `<table class="table table-dense">` +
                `<theader>${header}</theader>` +
                `<tbody>${html}<tbody>` +
                `</table>`;
        })
        .catch((e) => {
            console.log(`Got an error ${e}`);
        });
});

