<script lang="ts">
  import DataTableSearch from "$lib/components/DataTable/DataTableSearch.svelte";
  import Navbar from "$lib/components/Navbar.svelte";

  const { study }: StudyProps = $props();

  const datasetColumns = [
      { key: "name", label: "Dataset" },
      { key: "year", label: "Year" },
      { key: "region", label: "Region" },
      { key: "PMID", label: "PMID" },
      { key: "species", label: "Species" },
      { key: "author", label: "Author" },
      { key: "disease", label: "Disease" },
    ];
</script>

<svelte:head>
  <title>MultiSC-Viewer Study</title>
</svelte:head>

<Navbar />
<div style="margin:auto 20 rem;">
  <h1>{study.name}</h1>
  <p>{study.description}</p>
</div>

<DataTableSearch label="Datasets" data={study.datasets} columns={datasetColumns} />
