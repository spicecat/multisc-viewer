<script lang="ts">
  import Navbar from "$lib/components/Navbar.svelte";
  import DatasetSearch from "$lib/components/DatasetSearch.svelte";

  const { publication }: PublicationProps = $props();
</script>

<svelte:head>
  <title>MultiSC-Viewer Publication</title>
</svelte:head>

<Navbar />
<main style="margin:auto 10rem;">
  <section>
    <h1>{publication.name}</h1>
    <p>{publication.description}</p>
  </section>

  <hr />

  <section>
    <h1>Datasets</h1>
    <DatasetSearch datasets={publication.datasets} />
  </section>
</main>
