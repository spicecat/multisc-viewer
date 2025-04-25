<script lang="ts">
  import DataTableSearch from "$lib/components/DataTable/DataTableSearch.svelte";
  import Navbar from "$lib/components/Navbar.svelte";
  import { onMount } from "svelte";
  import type { Publication } from "../../interfaces/types";

  // State
  let publications: Publication[] = $state([]);

  const publicationColumns = [
    { key: "publication", label: "Publication", url: true },
    { key: "description", label: "Description" },
    { key: "year", label: "Year" },
    { key: "species", label: "Species" },
    { key: "author", label: "Authors" },
    { key: "disease", label: "Disease" },
    { key: "cellType", label: "Cell Type" },
    { key: "datasets", label: "Datasets" },
  ];

  // Fetch and process publication data
  onMount(() => {
    fetch("/publication")
      .then((res) => res.json())
      .then((data) => {
        publications = data.map((publication) => ({
          ...publication,
          publication: { name: publication.name, href: `/publication/${publication.publicationId}` },
          year: [
            ...new Set(publication.datasets.map((dataset) => dataset.year)),
          ].join(" "),
          species: [
            ...new Set(publication.datasets.map((dataset) => dataset.species)),
          ].join(" "),
          author: [
            ...new Set(publication.datasets.map((dataset) => dataset.author)),
          ].join(" "),
          disease: [
            ...new Set(publication.datasets.map((dataset) => dataset.disease).flat()),
          ].join(" "),
          cellType: [
            ...new Set(publication.datasets.map((dataset) => dataset.cellType)),
          ].join(" "),
          datasets: publication.datasets.length,
        }));
      });
  });
</script>

<svelte:head>
  <title>MultiSC-Viewer - Home</title>
</svelte:head>

<Navbar />

<main style="margin:auto 10rem;">
  <section>
    <h1>MultiSC-Viewer</h1>
    <p>
      MultiSC-Viewer is a web-based tool designed to compare and visualize gene
      expression in multiple single cell/nucleus dataset side by side for any
      gene across different brain regions, disease conditions, and species.
    </p>
  </section>

  <section>
    <h2>Key Features</h2>
    <p>
      Multiple datasets visualization and comparison: MultiSC-Viewer is the only
      tool current available that supports a side-by-side comparative view of
      more than two datasets simultaneously. Interactive Data Visualization:
      Utilize dynamic, interactive graphs powered by ??Plotly for a real-time
      exploration of datasets, facilitating a deeper understanding of cellular
      diversity and expression similarities and differences across datasets.
      This framework uniquely combines Svelte with Nest.js: MultiSC-Viewer
      Server-side rendering: Page requests are handled by the server, which
      pre-computes some page layout and data elements before responding,
      speeding up loading and reducing latency and allowing faster/more
      controlled access to certain types of data.
    </p>
  </section>

  <section>
    <h2>Local Installation</h2>
    <p>
      MultiSC-Viewer can be installed locally, supporting Linux, macOS, and
      Windows operating systems. This allows researchers to run the platform in
      their own environments, facilitating private data management and making
      use of local computational resources. For detailed installation
      instructions and access to additional resources, please visit Link.
    </p>
  </section>
  
  <hr />
  
  <section>
    <h1>Publications</h1>
    <DataTableSearch label="Publications" data={publications} columns={publicationColumns} />
  </section>
</main>
