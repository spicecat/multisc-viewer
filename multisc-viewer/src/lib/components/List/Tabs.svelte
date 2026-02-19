<script lang="ts" generics="T">
import { Tabs } from "@skeletonlabs/skeleton-svelte";
import type { Snippet } from "svelte";

const {
	tabs,
	data,
	filter,
	panels,
}: {
	tabs: string[];
	data: T[];
	filter: (tab: string) => (data: T) => boolean;
	panels: Snippet<[T[]]>;
} = $props();

let tab = $state("");

const items = $derived(tab ? data.filter(filter(tab)) : data);
</script>

<Tabs value={tab} onValueChange={(e) => (tab = e.value)}>
  <Tabs.List>
    <Tabs.Trigger class="flex-1" value="">
      <span class="font-bold">{data.length}</span>
      All
    </Tabs.Trigger>
    {#each tabs as tab (`tab-${tab}`)}
      <Tabs.Trigger class="flex-1" value={tab}>
        <span class="font-bold">{data.filter(filter(tab)).length}</span>
        {tab}
      </Tabs.Trigger>
    {/each}
  </Tabs.List>
  <Tabs.Content value={tab}>
    {@render panels(items)}
  </Tabs.Content>
</Tabs>
