<script lang="ts">
  import DataTableSearch from "$lib/components/DataTable/DataTableSearch.svelte";
  import Navbar from "$lib/components/Navbar.svelte";

  let studies = [
    {
      study: { url: "study1", name: "Example Study 1" },
      authors: "Author A, Author B",
      diseases: "Disease X, Disease Y",
    },
    {
      study: { url: "study2", name: "Example Study 2" },
      authors: "Author E",
      diseases: "Disease X, Disease Z",
    },
  ];

  const studyColumns = [
    { key: "study", label: "Study", url: true },
    { key: "authors", label: "Authors" },
    { key: "diseases", label: "Diseases" },
  ];
</script>

<svelte:head>
  <title>Dataset Comparison About</title>
</svelte:head>

<Navbar />
<div style="margin:1rem;">
  <h1>About</h1>
  <p>This is a website to compare plots across datasets.</p>
</div>

<DataTableSearch label="Studies" data={studies} columns={studyColumns} />
