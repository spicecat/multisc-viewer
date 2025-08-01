<script lang="ts">
	let { geneControlsOpen = $bindable() }: { geneControlsOpen: boolean } = $props();
</script>

<button
	class="btn-icon btn"
	onclick={() => (geneControlsOpen = !geneControlsOpen)}
	aria-label={geneControlsOpen ? 'Close gene controls' : 'Open gene controls'}
>
	<span class="material-icons">{geneControlsOpen ? 'close' : 'menu'}</span>
</button>
