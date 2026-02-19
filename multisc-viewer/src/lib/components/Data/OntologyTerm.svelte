<script lang="ts">
import { Globe } from "@lucide/svelte";
import type { Snippet } from "svelte";
import Chip from "$lib/components/List/Chip.svelte";
import type { OntologyTerm } from "$lib/types/daemon";

let {
	ontologyTerm = [],
	tags = $bindable(),
	tag,
	children,
}: {
	ontologyTerm?: OntologyTerm;
	tags?: string[];
	tag?: string;
	children?: Snippet;
} = $props();
</script>

<span class="inline-flex flex-wrap gap-1">
	{@render children?.()}
	{#each ontologyTerm as ot (`chip-ot-${ot.name}`)}
		{@const name = ot.displayName ?? ot.name}
		<Chip bind:tags tag={tag && `${tag}=${name}`}>
			{name}
			{#if ot.url}
				<!-- eslint-disable svelte/no-navigation-without-resolve -->
				<a
					class="anchor"
					href={ot.url}
					target="_blank"
					title={`${name} ontology source`}
					rel="external noopener noreferrer"
				>
					<Globe size="16" />
				</a>
				<!-- eslint-enable svelte/no-navigation-without-resolve -->
			{/if}
		</Chip>
	{/each}
</span>
