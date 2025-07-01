<script lang="ts">
  import Navbar from "$lib/components/Navbar.svelte";
  let { children } = $props();
</script>

<Navbar />

<main style="margin:auto 10rem;">
  {@render children()}
</main>
