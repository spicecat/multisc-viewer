<script lang="ts">
  import Button, { Label } from "@smui/button";

  let { geneControlsOpen = $bindable() } = $props();
</script>

<Button onclick={() => (geneControlsOpen = !geneControlsOpen)} variant="raised">
  <Label>Gene Controls</Label>
</Button>
