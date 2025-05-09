<script lang="ts">
  import DatasetsSelector from "$lib/components/DatasetsSelector.svelte";
  import Navbar from "$lib/components/Navbar.svelte";

  let { datasets }: { datasets: Dataset[] } = $props();
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
      <a href="https://git.jasonxu.dev/JasonXu/plot-viewer/">Source Code</a>
    </p>
  </section>

  <section>
    <h2>Key Features</h2>
    <p>
      Multiple datasets visualization and comparison: MultiSC-Viewer is the only
      tool current available that supports a side-by-side comparative view of
      more than two datasets simultaneously.
    </p>
  </section>

  <hr />

  <section>
    <h2>Datasets</h2>
    <DatasetsSelector {datasets} />
  </section>
</main>
