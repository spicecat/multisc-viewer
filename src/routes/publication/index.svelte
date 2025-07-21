<script lang="ts">
  import PublicationsTable from "$lib/components/DataTable/PublicationsTable.svelte";
  import Navbar from "$lib/components/Navbar.svelte";

  let { publications }: { publications: Publication[] } = $props();
</script>

<svelte:head>
  <title>MultiSC-Viewer - Publications</title>
</svelte:head>

<Navbar />

<main style="margin:auto 10rem;">
  <section>
    <h1>Publications</h1>
    <PublicationsTable {publications} />
  </section>
</main>
