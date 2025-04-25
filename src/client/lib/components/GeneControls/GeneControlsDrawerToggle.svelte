<script lang="ts">
  import IconButton from "@smui/icon-button";

  // Props
  let { geneControlsOpen = $bindable() } = $props();
</script>

<IconButton 
  class="material-icons toggle-button" 
  onclick={() => geneControlsOpen = !geneControlsOpen}
  aria-label={geneControlsOpen ? "Close gene controls" : "Open gene controls"}
>
  {geneControlsOpen ? 'close' : 'menu'}
</IconButton>

<style lang="scss">
  .toggle-button {
    margin-right: 0.5rem;
  }
</style>
