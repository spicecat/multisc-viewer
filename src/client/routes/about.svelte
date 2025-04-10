<script lang="ts">
  import DataTableSearch from "$lib/components/DataTable/DataTableSearch.svelte";
  import Navbar from "$lib/components/Navbar.svelte";
  import { onMount } from "svelte";

  let studies: Study[] = $state([]);

  const studyColumns = [
    { key: "study", label: "Study", url: true },
    { key: "description", label: "Description" },
  ];

  onMount(() => {
    fetch("/study")
      .then((res) => res.json())
      .then((data) => {studies = data.map(study=>
        ({...study, study: {name: study.name, href: `/study/${study.name}`}})
      );});
  });
</script>

<svelte:head>
  <title>MultiSC-Viewer About</title>
</svelte:head>

<Navbar />
<div style="margin:auto 20rem;">
  <h1>About</h1>
  <p>
    MultiSC-Viewer is a website to compare and visualize gene expression in
    multiple single cell/nucleus datasets across different brain regions,
    disease conditions, and species.
  </p>
  <hr />
  <h1>Studies</h1>
  <DataTableSearch label="Studies" data={studies} columns={studyColumns} />
</div>
