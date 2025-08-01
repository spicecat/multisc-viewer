<script lang="ts">
	import { RadioGroup, RadioItem } from '@skeletonlabs/skeleton-svelte';

	let { groupBy = $bindable() }: { groupBy: Grouping } = $props();

	const groupOptions = [
		{ value: 'Genotype', label: 'Group by Genotype' },
		{ value: 'CellType', label: 'Group by Cell Type' }
	];
</script>

<RadioGroup>
	{#each groupOptions as { value, label }}
		<RadioItem bind:group={groupBy} name="groupBy" {value}>{label}</RadioItem>
	{/each}
</RadioGroup>
