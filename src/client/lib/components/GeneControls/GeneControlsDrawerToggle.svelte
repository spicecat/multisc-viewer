<script lang="ts">
  import IconButton from "@smui/icon-button";

  let { geneControlsOpen = $bindable() }: {geneControlsOpen: boolean} = $props();
</script>

<IconButton
  class="material-icons"
  onclick={() => (geneControlsOpen = !geneControlsOpen)}
  aria-label={geneControlsOpen ? "Close gene controls" : "Open gene controls"}
>
  {geneControlsOpen ? "close" : "menu"}
</IconButton>
