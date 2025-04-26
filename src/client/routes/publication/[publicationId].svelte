<script lang="ts">
  import Navbar from "$lib/components/Navbar.svelte";
  import DatasetsSelector from "$lib/components/DatasetsSelector.svelte";

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
    <DatasetsSelector datasets={publication.datasets} />
  </section>
</main>
